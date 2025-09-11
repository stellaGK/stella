!###############################################################################
!                     INITIALISE THE DISTRIBUTION FUNCTION                      
!###############################################################################
! 
! The subroutine init_arrays_distribution_function() will allocate the following
! arrays, and initialise their values to zero:
!     - gnew(naky, nakx, -nzgrid:nzgrid, ntubes, -vmu-layout-)
!     - gold(naky, nakx, -nzgrid:nzgrid, ntubes, -vmu-layout-)
!     - g_scratch(naky, nakx, -nzgrid:nzgrid, ntubes, -vmu-layout-)
!     - gvmu(nvpa, nmu, -kxkyz-layout-)
! 
! If <split_parallel_dynamics> = True the following array will be initialised:
!     - g_kymus(nakx, -nzgrid:nzgrid, ntubes, vpa, -kymus-layout-)
! 
! The distribution function is initialised by initialising first the electrostatic
! potential <phi> according to one of the following options:
!     - Maxwellian distribution
!     - Noise based on a Random Number Generator associated to <rng_seed>
!     - Read the saved distribution function from an old simulation, so that
!       we can restart the simulation and continue the time advance
!     - kpar ...
!     - rh ...
!     - rmap ...
! 
! Once the electrostatic potential is initialised it is used to initialise 
! <gvmu> = g(mu,vpa,ikxkys) in the init_distribution_function_vs_muvpa() routine.
! 
! Next the distribution function on the ikxkys-layout is scattered to the ivpamus-
! layout to initialise <gold> and <gnew> defined as g(kx,ky,z,ivpamus) through the
! init_distribution_function_vs_kxkyz() routine.
! 
!###############################################################################
module initialise_distribution_function

   ! Load debug flags
   use debug_flags, only: debug => dist_fn_debug

   implicit none
   
   ! Public routines
   public :: read_parameters_distribution_function
   public :: init_distribution_function
   public :: finish_distribution_function
   public :: reset_init
   
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

   ! Choose the initalization option for the potential
   integer :: init_distribution_switch
   
   ! This variable is read when the module is initialised and used later in stella_restore()
   real :: scale
   
   ! During the initialization of this module we set the restart path
   ! and we will parse it to stella_save.fpp through save_init()
   character(len=300) :: restart_file
   character(len=150) :: restart_dir

   ! Remember whether the routines have already been initialised
   logical :: initialised_read_parameters = .false.
   logical :: initialised_arrays = .false.
   logical :: initialised_distribution_function = .false.
   logical :: initialised_distribution_function_vs_muvpa = .false.
   logical :: initialised_distribution_function_vs_kxkyz = .false.

contains


!###############################################################################
!############################### READ INPUT FILE ###############################
!###############################################################################

   !****************************************************************************
   !                              READ INPUT FILE                              !
   !****************************************************************************
   subroutine read_parameters_distribution_function

      use mp, only: proc0, broadcast
      use stella_save, only: init_save, read_many
      use stella_layouts, only: read_parameters_parallelisation_layouts
      use system_fortran, only: systemf
      use stella_save, only: read_many
      
      ! Read namelist from input file
      use namelist_initialise_distribution_function, only: read_namelist_initialise_distribution
      use namelist_initialise_distribution_function, only: read_namelist_restart_options
      use namelist_initialise_distribution_function, only: read_namelist_initialise_distribution_noise

      ! Load the <init_distribution_switch> parameters
      use namelist_initialise_distribution_function, only: init_distribution_option_maxwellian
      use namelist_initialise_distribution_function, only: init_distribution_option_noise
      use namelist_initialise_distribution_function, only: init_distribution_option_restart_many
      use namelist_initialise_distribution_function, only: init_distribution_option_kpar
      use namelist_initialise_distribution_function, only: init_distribution_option_rh
      use namelist_initialise_distribution_function, only: init_distribution_option_remap

      implicit none

      integer :: ind_slash
      
      !-------------------------------------------------------------------------

      ! Only initialise once
      if (initialised_read_parameters) return
      initialised_read_parameters = .true.
   
      ! Make sure the parallelisation parameters are already read
      call read_parameters_parallelisation_layouts
      
      ! Read <initialise_distribution> namelist
      call read_namelist_initialise_distribution(init_distribution_switch, phiinit, scale_to_phiinit)
         
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
      if (ind_slash == 0) then ! No slash present
         restart_file = trim(restart_dir)//trim(restart_file)
      else ! Slash present
         restart_file = trim(restart_file(1:ind_slash))//trim(restart_dir)//trim(restart_file(ind_slash + 1:))
      end if 
      
      ! Initialize the netcdf saving
      call init_save(restart_file)

   end subroutine read_parameters_distribution_function
   
!###############################################################################
!##################### INITIALISE THE DISTRIBUTION FUNCTION ####################
!###############################################################################

   !****************************************************************************
   !                   INITIALISE THE DISTRIBUTION FUNCTION                    !
   !****************************************************************************
   subroutine init_distribution_function(restarted, istep0)
      
      implicit none
      
      ! Arguments
      logical, intent(out) :: restarted
      integer, intent(out) :: istep0
      
      !-------------------------------------------------------------------------
      
      ! Only initialise once
      if (initialised_distribution_function) return
      initialised_distribution_function = .true.
      
      ! Allocate the distribution-sized arrays and initialise them to zero
      if (debug) write (6, *) "stella::init_stella::init_arrays_distribution_function"
      call init_arrays_distribution_function
      
      ! Initialise the guiding-center distribution function <gvmu>(nvpa, nmu, -kxkyzs-layout-)
      if (debug) write (6, *) "stella::init_stella::init_distribution_function_vs_muvpa"
      call init_distribution_function_vs_muvpa(restarted, istep0)
      
      ! Initialise the guiding-center distribution function <gnew>(kx, ky, z, -vpamus-layout-)
      ! Use mapping from kxkyz_lo to vmu_lo to get a copy of g that has ky, kx and z local to each core;
      ! This distribution function is stored in <gnew> and copied to <gold>
      if (debug) write (6, *) "stella::init_stella::init_distribution_function_vs_kxkyz"
      call init_distribution_function_vs_kxkyz(restarted)
      
   end subroutine init_distribution_function

   !****************************************************************************
   !                              ALLOCATE ARRAYS                               
   !****************************************************************************
   ! This subroutine initialises the distribution function arrays used in the
   ! STELLA code. It allocates the arrays and sets them to zero.
   !****************************************************************************
   subroutine init_arrays_distribution_function

      ! Parallelisation
      use stella_layouts, only: kxkyz_lo
      use stella_layouts, only: vmu_lo
      use stella_layouts, only: kymus_lo
      
      ! Distribution functions that are allocated here
      use arrays_distribution_function, only: gnew, gold, g_scratch
      use arrays_distribution_function, only: gvmu, g_kymus
      use arrays_distribution_function, only: g0, g1, g2, g3
      
      ! Numerical flags
      use parameters_numerical, only: split_parallel_dynamics
      use parameters_numerical, only: explicit_algorithm_switch
      use parameters_numerical, only: explicit_algorithm_rk4
      
      ! Grids
      use grids_z, only: nzgrid, ntubes
      use grids_kxky, only: naky, nakx
      use grids_velocity, only: nvpa, nmu
      
      implicit none
      
      !-------------------------------------------------------------------------

      ! Only initialise once
      if (initialised_arrays) return
      initialised_arrays = .true.

      ! Allocate arrays
      if (debug) write (*, *) 'dist_fn::init_arrays_distribution_function::allocate_arrays'
      if (.not. allocated(gnew)) allocate (gnew(naky, nakx, -nzgrid:nzgrid, ntubes, vmu_lo%llim_proc:vmu_lo%ulim_alloc))
      if (.not. allocated(gold)) allocate (gold(naky, nakx, -nzgrid:nzgrid, ntubes, vmu_lo%llim_proc:vmu_lo%ulim_alloc))
      if (.not. allocated(g_scratch)) allocate (g_scratch(naky, nakx, -nzgrid:nzgrid, ntubes, vmu_lo%llim_proc:vmu_lo%ulim_alloc))
      if (.not. allocated(gvmu)) allocate (gvmu(nvpa, nmu, kxkyz_lo%llim_proc:kxkyz_lo%ulim_alloc))
         
      ! Initialise arrays to zero
      gnew = 0.
      gold = 0.
      g_scratch = 0.
      gvmu = 0.
      
      !-------------------------------------------------------------------------
      
      ! Only allocate <g_kymus> if <split_parallel_dynamics> = True
      if (.not. allocated(g_kymus)) then
            if (.not. split_parallel_dynamics) then
               allocate (g_kymus(nakx, -nzgrid:nzgrid, ntubes, nvpa, kymus_lo%llim_proc:kymus_lo%ulim_alloc))
            else
               allocate (g_kymus(1, 1, 1, 1, 1))
            end if
            g_kymus = 0.
      end if
      
      !-------------------------------------------------------------------------

      ! Initialise dummy arrays used inside the Runge-Kutta schemes and diagnostics
      if (.not. allocated(g0)) allocate (g0(naky, nakx, -nzgrid:nzgrid, ntubes, vmu_lo%llim_proc:vmu_lo%ulim_alloc))
      if (.not. allocated(g1)) allocate (g1(naky, nakx, -nzgrid:nzgrid, ntubes, vmu_lo%llim_proc:vmu_lo%ulim_alloc))
      if (.not. allocated(g2)) allocate (g2(naky, nakx, -nzgrid:nzgrid, ntubes, vmu_lo%llim_proc:vmu_lo%ulim_alloc))
      
      ! Initialise arrays to zero
      g0 = 0.
      g1 = 0.
      g2 = 0.
      
      ! Only allocate <g3> if the 4th order Runge-Kutta scheme is utilised
      if (.not. allocated(g3)) then
         if (explicit_algorithm_switch == explicit_algorithm_rk4) then
            allocate (g3(naky, nakx, -nzgrid:nzgrid, ntubes, vmu_lo%llim_proc:vmu_lo%ulim_alloc))
         else
            allocate (g3(1, 1, 1, 1, 1))
         end if
         g3 = 0.
      end if
      
   end subroutine init_arrays_distribution_function
   
   !****************************************************************************
   !       INITIALISE THE DISTRIBUTION FUNCTION <GVMU> = G(MU,VPA,IKXKYZ)      !
   !****************************************************************************
   subroutine init_distribution_function_vs_muvpa(restarted, istep0)

      ! Flags
      use stella_save, only: init_tstart
      use parameters_numerical, only: maxwellian_normalization
      
      ! Load the <init_distribution_switch> parameters
      use namelist_initialise_distribution_function, only: init_distribution_option_maxwellian
      use namelist_initialise_distribution_function, only: init_distribution_option_noise
      use namelist_initialise_distribution_function, only: init_distribution_option_restart_many
      use namelist_initialise_distribution_function, only: init_distribution_option_kpar
      use namelist_initialise_distribution_function, only: init_distribution_option_rh
      use namelist_initialise_distribution_function, only: init_distribution_option_remap

      ! Arguments
      logical, intent(out) :: restarted
      integer, intent(out) :: istep0
      
      ! Local variables
      integer :: istatus

      !-------------------------------------------------------------------------

      ! Only initialise once
      if (initialised_distribution_function_vs_muvpa) return
      initialised_distribution_function_vs_muvpa = .false.

      ! Assume this is a new simulation, starting from time step <istep0> = 0.
      ! If <init_distribution_switch> = <init_distribution_option_restart_many>,
      ! then we have restarted the simulation from <istep0> > 0.
      restarted = .false.
      istep0 = 0
      
      ! Initialise the distribution function <gvmu> = g(mu,vpa,ikxkys) based on 
      ! one of the following options, which has been selected in the input file
      select case (init_distribution_switch)
      case (init_distribution_option_maxwellian)
         call initialise_distribution_maxwellian
      case (init_distribution_option_noise)
         call initialise_distribution_noise
      case (init_distribution_option_kpar)
         call initialise_distribution_kpar
      case (init_distribution_option_rh)
         call initialise_distribution_rh
      case (init_distribution_option_remap)
         call initialise_distribution_remap
      case (init_distribution_option_restart_many)
         call initialise_distribution_restart_many
         call init_tstart(tstart, istep0, istatus)
         restarted = .true.
         scale = 1.
      end select

      ! If <maxwwellian_normalization> = .true., the pdf is normalized by F0 (which is not the case otherwise)
      ! unless reading in g from a restart file, normalise g by F0 for a full flux surface simulation
      if (maxwellian_normalization .and. init_distribution_switch /= init_distribution_option_restart_many) then
         call normalize_by_maxwellian
      end if

   end subroutine init_distribution_function_vs_muvpa

   !*****************************************************************************
   ! INITIALISE THE DISTRIBUTION FUNCTIONS <GOLD> = <GNEW> = G(KX,KY,Z,IVPAMUS) !
   !*****************************************************************************
   ! This subroutine initialises the <gold> and <gnew> = g(kx,ky,z,ivpamus)
   ! arrays, based on the previously initialised <gvmu> = g(mu,vpa,ikxkys) array.
   !****************************************************************************
   subroutine init_distribution_function_vs_kxkyz(restarted)
   
      ! Parallelisation
      use redistribute, only: gather, scatter
      use calculations_redistribute, only: kxkyz2vmu

      ! Distribution function
      use arrays_distribution_function, only: gvmu, gold, gnew
      
      ! Physics flags
      use parameters_physics, only: radial_variation

      implicit none

      ! Arguments
      logical, intent(in) :: restarted
      
      !-------------------------------------------------------------------------

      ! Only initialise once
      if (initialised_distribution_function_vs_kxkyz) return
      initialised_distribution_function_vs_kxkyz = .false.

      ! We have already initialised <gvmu> = g(mu,vpa,ikxkys) in the init_distribution_
      ! function_vs_muvpa() routine. Now we get <gnew> = g(kx,ky,z,ivpamus).
      call gather(kxkyz2vmu, gvmu, gnew)

      ! Calculate radial corrections to F0 for use in the Krook operator.
      ! Note that this routine will redefine <gnew>.
      if (radial_variation) call add_corrections_to_g_for_radial_variation(restarted)

      ! Initialise <gold> to be a copy of <gnew>
      gold = gnew
      
   end subroutine init_distribution_function_vs_kxkyz
   
!###############################################################################
!############################## RADIAL VARIATION ###############################
!###############################################################################
! Calculate the radial corrections to F0 for use in the Krook operator.
! It also initialises the gold array with the values of gnew.
! If the radial variation is not enabled, it simply copies gnew to gold.
! The gxyz arrays are used in the Krook operator and projection method.
!###############################################################################
    
   subroutine add_corrections_to_g_for_radial_variation(restarted)

      ! Parallelisation
      use redistribute, only: scatter
      use calculations_redistribute, only: kxkyz2vmu
      use stella_layouts, only: vmu_lo, iv_idx, imu_idx, is_idx
      
      ! Distribution function
      use arrays_distribution_function, only: gvmu, gnew
   
      ! Calculations
      use calculations_transforms, only: transform_kx2x_xfirst, transform_x2kx_xfirst
      use calculations_kxky, only: multiply_by_rho
      
      ! Grids
      use grids_kxky, only: nalpha, nakx, naky
      use grids_velocity, only: mu, vpa, vperp2
      use grids_z, only: nzgrid, ntubes
      use grids_species, only: spec, pfac
      
      ! Geometry
      use geometry, only: dBdrho, gfac

      implicit none

      ! Arguments
      logical, intent(in) :: restarted
   
      ! Local variables
      real :: corr
      integer :: ivmu, is, imu, iv, it, iz, ia
      real, dimension(:, :), allocatable :: energy
      complex, dimension(:, :), allocatable :: g0k
      
      !----------------------------------------------------------------------

      ! Assume we only have a single field line
      ia = 1
   
      ! Allocate local arrays
      allocate (energy(nalpha, -nzgrid:nzgrid))
      allocate (g0k(naky, nakx))

      ! Iterate over (mu,vpa,s)
      do ivmu = vmu_lo%llim_proc, vmu_lo%ulim_proc
         is = is_idx(vmu_lo, ivmu)
         imu = imu_idx(vmu_lo, ivmu)
         iv = iv_idx(vmu_lo, ivmu)
         energy = (vpa(iv)**2 + vperp2(:, :, imu)) * (spec(is)%temp_psi0 / spec(is)%temp)
         
         do it = 1, ntubes
            do iz = -nzgrid, nzgrid

               corr = -(pfac * (spec(is)%fprim + spec(is)%tprim * (energy(ia, iz) - 1.5)) &
                        + 2 * gfac * mu(imu) * dBdrho(iz))

               if (.not. restarted) then
                  g0k = corr * gnew(:, :, iz, it, ivmu)
                  call multiply_by_rho(g0k)
                  gnew(:, :, iz, it, ivmu) = gnew(:, :, iz, it, ivmu) + g0k
               end if
               
            end do
         end do
      end do
      
      ! Deallocate local arrays
      deallocate (energy, g0k)

      ! Recalculate <gvmu> based on the new <gnew>
      if (.not. restarted) call scatter(kxkyz2vmu, gnew, gvmu)
   
   end subroutine add_corrections_to_g_for_radial_variation
      
!###############################################################################
!############### OPTIONS TO INITIALISE <GVMU> = G(MU,VPA,IKXKYZ) ###############
!###############################################################################

   !****************************************************************************
   !                     INITIALISE POTENTIAL: MAXWELLIAN                      !
   !****************************************************************************
   subroutine initialise_distribution_maxwellian

      use mp, only: proc0, broadcast
      use constants, only: zi
      use grids_species, only: spec
      use grids_z, only: nzgrid, zed
      use grids_kxky, only: naky, nakx, ikx_max
      use grids_kxky, only: theta0, akx, zonal_mode
      use grids_kxky, only: reality
      use grids_velocity, only: nvpa, nmu
      use grids_velocity, only: vpa
      use grids_velocity, only: maxwell_vpa, maxwell_mu, maxwell_fac
      use arrays_distribution_function, only: gvmu
      use stella_layouts, only: kxkyz_lo, iz_idx, ikx_idx, iky_idx, is_idx
      use ran, only: ranf
      use namelist_initialise_distribution_function, only: read_namelist_initialise_distribution_maxwellian

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

   end subroutine initialise_distribution_maxwellian

   !****************************************************************************
   !                       INITIALISE POTENTIAL: NOISE                         !
   !****************************************************************************
   subroutine initialise_distribution_noise

      use mp, only: proc0, broadcast
      use arrays, only: kperp2
      use grids_species, only: spec
      use grids_z, only: nzgrid, ntubes
      use grids_extended_zgrid, only: ikxmod, nsegments, neigen
      use grids_extended_zgrid, only: it_right
      use grids_extended_zgrid, only: periodic, phase_shift
      use grids_kxky, only: naky, nakx, reality
      use grids_kxky, only: zonal_mode
      use grids_velocity, only: nvpa, nmu
      use grids_velocity, only: maxwell_vpa, maxwell_mu, maxwell_fac
      use arrays_distribution_function, only: gvmu
      use stella_layouts, only: kxkyz_lo
      use stella_layouts, only: iky_idx, ikx_idx, iz_idx, it_idx, is_idx
      use mp, only: proc0, broadcast, max_allreduce
      use mp, only: scope, crossdomprocs, subprocs
      use file_utils, only: runtype_option_switch, runtype_multibox
      use parameters_physics, only: include_nonlinear
      use ran
      use namelist_initialise_distribution_function, only: read_namelist_initialise_distribution_noise

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
         call initialise_distribution_maxwellian
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
      end if

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

   end subroutine initialise_distribution_noise

   !****************************************************************************
   !                       INITIALISE POTENTIAL: KPAR                          !
   !****************************************************************************
   subroutine initialise_distribution_kpar 
   
      use mp, only: proc0, broadcast
      use grids_z, only: nzgrid, zed
      use grids_kxky, only: naky, nakx
      use grids_kxky, only: theta0
      use grids_velocity, only: nvpa, nmu
      use grids_velocity, only: vpa, vperp2
      use grids_velocity, only: maxwell_vpa, maxwell_mu, maxwell_fac
      use arrays_distribution_function, only: gvmu
      use stella_layouts, only: kxkyz_lo, iky_idx, ikx_idx, iz_idx, is_idx
      use namelist_initialise_distribution_function, only: read_namelist_initialise_distribution_kpar
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

   end subroutine initialise_distribution_kpar

   !****************************************************************************
   !                        INITIALISE POTENTIAL: RH                           !
   !****************************************************************************
   subroutine initialise_distribution_rh

      use mp, only: proc0, broadcast
      use grids_species, only: spec
      use arrays_distribution_function, only: gvmu
      use arrays, only: kperp2
      use stella_layouts, only: kxkyz_lo
      use stella_layouts, only: iky_idx, ikx_idx, iz_idx, is_idx
      use grids_velocity, only: maxwell_vpa, maxwell_mu, maxwell_fac
      use grids_velocity, only: nvpa, nmu
      use grids_kxky, only: akx
      use namelist_initialise_distribution_function, only: read_namelist_initialise_distribution_rh

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

      ! initialise g to be a Maxwellian with a constant density perturbation

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

   end subroutine initialise_distribution_rh

   !****************************************************************************
   !                      INITIALISE POTENTIAL: REMAP                          !
   !****************************************************************************
   subroutine initialise_distribution_remap

      use grids_species, only: spec
      use arrays_distribution_function, only: gvmu
      use stella_layouts, only: kxkyz_lo
      use stella_layouts, only: iky_idx, ikx_idx, iz_idx, is_idx
      use grids_velocity, only: maxwell_vpa, maxwell_mu, maxwell_fac
      use grids_velocity, only: nvpa, nmu

      implicit none

      integer :: ikxkyz, iky, ikx, iz, is, ia

      !-------------------------------------------------------------------------

      ! initialise g to be a Maxwellian with a constant density perturbation

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

   end subroutine initialise_distribution_remap

   !****************************************************************************
   !                       INITIALISE POTENTIAL: MANY                          !
   !****************************************************************************
   subroutine initialise_distribution_restart_many

      use arrays_distribution_function, only: gvmu
      use stella_save, only: stella_restore
      use mp, only: proc0
      use file_units, only: unit_error_file
      
      implicit none

      integer :: istatus

      !-------------------------------------------------------------------------

      ! should really check if profile_variation=T here but need
      ! to move profile_variation to module that is accessible here
      call stella_restore(gvmu, scale, istatus)

      if (istatus /= 0) then
         if (proc0) write (unit_error_file, *) "Error reading file: ", trim(restart_file)
         gvmu = 0.
      end if

   end subroutine initialise_distribution_restart_many
  

!###############################################################################
!################################ CALCULATIONS #################################
!###############################################################################

   !****************************************************************************
   !                        NORMALIZE BY MAXWELLIAN                            !
   !****************************************************************************
   subroutine normalize_by_maxwellian

      use stella_layouts, only: kxkyz_lo, is_idx, iz_idx
      use arrays_distribution_function, only: gvmu
      use grids_velocity, only: nmu
      use grids_velocity, only: maxwell_mu, maxwell_vpa, maxwell_fac

      implicit none

      integer :: ia, imu
      integer :: ikxkyz, iz, is

      !-------------------------------------------------------------------------

      ! gvmu is initialised with a Maxwellian weighting for flux tube simulations,
      ! with the Maxwellian evaluated at ia = 1
      ! we are undoing that weighting here, so also need to use ia = 1
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
   
      use namelist_initialise_distribution_function, only: init_distribution_option_restart_many
      
      implicit none

      init_distribution_switch = init_distribution_option_restart_many

   end subroutine reset_init

!###############################################################################
!######################## FINISH DISTRIBUTION FUNCTION #########################
!###############################################################################

   subroutine finish_distribution_function

      use stella_save, only: finish_save
      use arrays_distribution_function, only: gnew, gold, gvmu
      use arrays_distribution_function, only: g_scratch, g_kymus
      use arrays_distribution_function, only: g0, g1, g2, g3

      implicit none

      ! Deallocate arrays
      if (allocated(gnew)) deallocate (gnew)
      if (allocated(gold)) deallocate (gold)
      if (allocated(g_scratch)) deallocate (g_scratch)
      if (allocated(gvmu)) deallocate (gvmu)
      if (allocated(g_kymus)) deallocate (g_kymus)
      if (allocated(g0)) deallocate (g0)
      if (allocated(g1)) deallocate (g1)
      if (allocated(g2)) deallocate (g2)
      if (allocated(g3)) deallocate (g3)

      ! Reset flags
      initialised_read_parameters = .false.
      initialised_arrays = .false.
      initialised_distribution_function_vs_muvpa = .false.
      initialised_distribution_function_vs_kxkyz = .false.
   
      call finish_save

   end subroutine finish_distribution_function

end module initialise_distribution_function
