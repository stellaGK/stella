!###############################################################################
!                   STELLA LAYOUTS FOR PARALLELISATION
!###############################################################################
! This modules defines the layout/division of the (x-y-z); (kx-y-z); (kx-ky-z)
! and (v-mu-s) grid on the processors, that are needed for the distribution functions.
! The layouts depend on the grid divisions along (x,y,kx,ky,z,mu,v,species,tubes).
! 
! First, read the parameters in the namelist <parallelisation> to determine
! the layout of the xyzs grid and the vms grid and broadcast it.
! 
! Divide the total number of points for each grid over the processors.
!       kxkyz_layout    (naky*nakx*nzed*ntubes*nspec - 1) points
!       kxyz_layout     (ny*nakx*nzed*ntubes*nspec - 1) points
!       xyz_layout      (ny*nx*nzed*ntubes*nspec - 1) points
!       vmu_layout      (2*nvgrid*nmu*nspec - 1) points
! 
! Each processor will be in charge of <blocksize> points ranging from <llim_proc>
! up to <ulim_proc>. The "is_idx", "ix_idx", "iy_idx", "iz_idx", "it_idx", ...
! routines will return the values of (x, y, kx, ky, z, mu, v, species, tubes)
! associated with the point <idx> which is <ikxkyz>; <ikxyz>; <ixyz> or <ivmu>.
! 
! The KXYZ Layout options are "xyzs", "xzys", "yxzs", "zxys", "yzxs", which
! decide the order of (kx, y, z, species, tubes) of the points <idx>. The VMU
! layout options are "vms" and "mvs" which decide the order of (species, mu, vgrid)
!       x --> [1, ..., nx or nakx]
!       y --> [1, ..., ny or naky]
!       z --> [1, ..., nzed = 2*nzgrid+1]
!       s --> [1, ..., nspec]
!       t --> [1, ..., ntubes]
!       v --> [1, ..., 2*nvgrid = nvpa]
!       m --> [1, ..., nmu]
!       s --> [1, ..., nspec]
! 
! The "proc_id" routine returns the id of the processor.
!       proc_id --> [0, ..., nproc-1]
! 
! The "idx_kxkyz"; "idx_kxyz"; "idx_xyz" or "idx_vmu" routine returns the id of the point.
!       idx_kxkyz    -->    point (iky,ikx,iz,it,is)
!       idx_kxyz     -->    point (iy,ikx,iz,it,is)
!       idx_xyz      -->    point (iy,ix,iz,it,is)
!       idx_vmu      -->    point (iv,imu,is)
! 
! The "idx_local" routine returns true if the point <idx> lies on the current
! processor <iproc> and returns false if the point <idx> does not lie on it.
! 
! The routines "kxkyzidx2vmuidx", "kxyzidx2vmuidx" and "xyzidx2vmuidx" allows us
! to find the correspondence between the points on the different layout types:
!        <ikxkyz> = point (iky,ikx,iz,it,is)
!        <ikxyz>  = point (iy, ikx,iz,it,is)
!        <ixyz>   = point (iy, ix, iz,it,is)
!        <ivmu>   = point (iv,imu,is)
! As well as the corresponding values of (iv,imu,is,iy,ix,iky,ikx,iz,it).
! 
!-------------------------------- Integer limit ---------------------------------
! 
! When running very high resolutions runs (high in kx and ky) the ikxkyz index 
! might reach the integer limit, which is 2,147,483,647 on 32 bit systems.
! We can define 64-bit integers which reach up to 9,223,372,036,854,775,807.
! This might have to be changed in stella in the future, use:
! integer, parameter :: huge = selected_int_kind(15)
! 
!--------------------------------- Input file ----------------------------------
! 
!&parallelisation
!   xyzs_layout = 'yxzs'
!   vms_layout = 'vms'
!   mat_gen = .false.
!   mat_read = .false.
!   lu_option = 'default'
!   fields_kxkyz = .false.
!/
! 
!###############################################################################
module parallelisation_layouts

   ! Debug flags
   use debug_flags, only: debug => parallelisation_debug

   ! Each parallelisation layout contains the following variables:
   ! - The dimensions (nzgrid, nzed, ntubes, nalpha, ny, naky, nx, nakx, nvgrid, nvpa, nmu, nspec)
   ! - The proccesor information (iproc, llim_world, ulim_world, llim_proc, ulim_proc, ulim_alloc, blocksize)
   use common_types, only: vmu_layout_type
   use common_types, only: kxkyz_layout_type
   use common_types, only: kxyz_layout_type
   use common_types, only: xyz_layout_type
   use common_types, only: kymus_layout_type
   
   implicit none

   ! Make routines accesible to other modules
   public :: finish_layouts
   public :: init_dist_fn_layouts
   public :: read_parameters_parallelisation_layouts
   
   !------------------------------ Input variables -----------------------------
   
   ! When parallelising real space or velocity space, we can choose the order of the indices through
   ! xyzs_layout = [yxzs, yzxs, xyzs, zyxs, zxys, xzys] with the default being 'yxzs'
   ! vms_layout = [vms, mvs] with the default being 'vms'
   public :: xyzs_layout
   public :: vms_layout
   
   ! By default, stella uses g(ky, kx, z, -vpamus-layout-) to calculate the fields
   ! However, the fields can be calculated using g(vpa, mu, -kxkysz-layout-) if <fields_kxkyz> = True
   public :: fields_kxkyz
   
   ! Flags to read and write the response matrix to a NetCDF file
   public :: mat_gen
   public :: mat_read
   
   ! Construct the response matrix on 1 processor (none),
   ! on all processors of 1 node (local) or across all nodes (global)
   public :: lu_option_switch
   public :: lu_option_local
   public :: lu_option_none
   public :: lu_option_global
   
   ! Print extra info with respect to the parallelisation
   public :: verbose
   
   !-------------------------- Parallelisation Layouts -------------------------
   
   ! The (kx,ky,z,tube) dimensions are local and the (vpa, mu, s) dimensions are parallelised.
   ! Therefore, (2*nvgrid*nmu*nspec - 1) grid points are divided over the processors.
   public :: vmu_lo
   
   ! The (mu, vpa) dimensions are local and the (kx, ky, z, tube, s) dimensions are parallelised.
   ! Therefore, (naky*nakx*nzed*ntubes*nspec - 1) grid points are divided over the processors.
   public :: kxkyz_lo
   
   ! The (mu, vpa) dimensions are local and the (kx, y, z, tube, s) dimensions are parallelised.
   ! Therefore, (ny*nakx*nzed*ntubes*nspec - 1) grid points are divided over the processors.
   public :: kxyz_lo
   
   ! The (mu, vpa) dimensions are local and the (x, y, z, tube, s) dimensions are parallelised.
   ! Therefore, (ny*nx*nzed*ntubes*nspec - 1) grid points are divided over the processors.
   public :: xyz_lo          ! Only used if include_parallel_nonlinearity = True
   
   ! This layout has the kx, z and vpa dimensions local on each processor,
   ! with ky, mu and species available to spread over processors
   public :: kymus_lo        ! Only used if split_parallel_dynamics = True

   ! Find the equivalent points between different parallelisation layouts
   public :: kxkyzidx2vmuidx
   public :: kxyzidx2vmuidx
   public :: xyzidx2vmuidx
   public :: kymusidx2vmuidx

   ! Find the index of (kx,ky,z,s,mu,vpa) of <ikxkyzs>, <ivpamus>, ...
   public :: ikx_idx, iky_idx, iz_idx, is_idx, imu_idx, iv_idx
   public :: idx_local, iy_idx, it_idx, idx
   
   ! The "proc_id" routine returns the id of the processor, with proc_id --> [0, ..., nproc-1]
   public :: proc_id

   private
   
   !----------------------------------------------------------------------------
   
   ! When parallelising real space or velocity space, we can choose the order of the indices through
   ! xyzs_layout = [yxzs, yzxs, xyzs, zyxs, zxys, xzys] with the default being 'yxzs'
   ! vms_layout = [vms, mvs] with the default being 'vms'
   character(len=4) :: xyzs_layout
   character(len=3) :: vms_layout
   character(len=5) :: kymus_layout

   ! Define the 5 possible parallelisation Layouts
   type(kxkyz_layout_type) :: kxkyz_lo
   type(kxyz_layout_type) :: kxyz_lo
   type(xyz_layout_type) :: xyz_lo
   type(vmu_layout_type) :: vmu_lo
   type(kymus_layout_type) :: kymus_lo

   ! Find the index of the <tube> on <ikxkyzs>, <ivpamus>, ...
   interface it_idx
      module procedure it_idx_kxkyz
      module procedure it_idx_kxyz
      module procedure it_idx_xyz
   end interface

   ! Find the index of <z> on <ikxkyzs>, <ivpamus>, ...
   interface iz_idx
      module procedure iz_idx_kxkyz
      module procedure iz_idx_kxyz
      module procedure iz_idx_xyz
   end interface

   ! Find the index of <vpa> on <ikxkyzs>, <ivpamus>, ...
   interface iv_idx
      module procedure iv_idx_vmu
   end interface

   ! Find the index of <ky> on <ikxkyzs>, <ivpamus>, ...
   interface iky_idx
      module procedure iky_idx_kxkyz
      module procedure iky_idx_kymus
   end interface

   ! Find the index of <y> on <ikxkyzs>, <ivpamus>, ...
   interface iy_idx
      module procedure iy_idx_kxyz
      module procedure iy_idx_xyz
   end interface

   ! Find the index of <kx> on <ikxkyzs>, <ivpamus>, ...
   interface ikx_idx
      module procedure ikx_idx_kxkyz
      module procedure ikx_idx_kxyz
   end interface

   ! Find the index of <x> on <ikxkyzs>, <ivpamus>, ...
   interface ix_idx
      module procedure ix_idx_xyz
   end interface

   ! Find the index of <mu> on <ikxkyzs>, <ivpamus>, ...
   interface imu_idx
      module procedure imu_idx_vmu
      module procedure imu_idx_kymus
   end interface

   ! Find the index of <s> on <ikxkyzs>, <ivpamus>, ...
   interface is_idx
      module procedure is_idx_kxkyz
      module procedure is_idx_kxyz
      module procedure is_idx_xyz
      module procedure is_idx_vmu
      module procedure is_idx_kymus
   end interface

   ! Find the processor id on <ikxkyzs>, <ivpamus>, ...
   interface proc_id
      module procedure proc_id_kxkyz
      module procedure proc_id_kxyz
      module procedure proc_id_xyz
      module procedure proc_id_vmu
      module procedure proc_id_kymus
   end interface

   ! Find the ID of the point on <ikxkyzs>, <ivpamus>, ...
   interface idx
      module procedure idx_kxkyz
      module procedure idx_kxyz
      module procedure idx_xyz
      module procedure idx_vmu
      module procedure idx_kymus
   end interface

   ! Find the ID of the point on <ikxkyzs>, <ivpamus>, ...
   interface idx_local
      module procedure idx_local_kxkyz, iz_local_kxkyz
      module procedure idx_local_kxyz, iz_local_kxyz
      module procedure idx_local_xyz, iz_local_xyz
      module procedure idx_local_vmu, iz_local_vmu
      module procedure idx_local_kymus, iz_local_kymus
   end interface 
   
   !------------------------------ Input variables -----------------------------

   ! By default, stella uses g(ky, kx, z, -vpamus-layout-) to calculate the fields
   ! However, the fields can be calculated using g(vpa, mu, -kxkysz-layout-) if <fields_kxkyz> = True
   logical :: fields_kxkyz
   
   ! Flags to read and write the response matrix to a NetCDF file
   logical :: mat_gen
   logical :: mat_read
   
   ! Print extra info to command prompt
   ! This variable is not broadcasted so it can only be true on proc0
   logical :: verbose = .false.
   
   ! Construct the response matrix on 1 processor (none),
   ! on all processors of 1 node (local) or across all nodes (global)
   integer :: lu_option_switch
   integer, parameter :: lu_option_none = 1
   integer, parameter :: lu_option_local = 2
   integer, parameter :: lu_option_global = 3

contains

!###############################################################################
!                                INITIALISATION                                 
!###############################################################################

   !****************************************************************************
   !                       INITIALISE THE STELLA LAYOUTS                        
   !****************************************************************************
   ! Read the parameters in the namelist <layouts_knobs> to determine
   ! the layout of the xyzs grid and the vms grid and broadcast it.
   !****************************************************************************
   subroutine read_parameters_parallelisation_layouts

      use mp, only: proc0
      !use mp, only: broadcast
      
      implicit none
      
      logical, save :: initialised = .false.
      
      !-------------------------------------------------------------------------

      ! Only initialise once
      if (initialised) return
      initialised = .true.

      ! Read input parameters and broadcast to all processors
      if (proc0) call read_parameters
      call broadcast_parameters

   end subroutine read_parameters_parallelisation_layouts

   !****************************************************************************
   !                               READ PARAMETERS                              
   !****************************************************************************
   ! Read the parameters in the namelist <parallelisation> to determine
   ! the layout of the xyzs grid and the vms grid.
   !****************************************************************************
   subroutine read_parameters

      use namelist_parallelisation, only: read_namelist_parallelisation
      
      implicit none

      call read_namelist_parallelisation(xyzs_layout, vms_layout, kymus_layout, &
          mat_gen, mat_read, lu_option_switch, fields_kxkyz, verbose)

   end subroutine read_parameters

   !****************************************************************************
   !                            BROADCAST PARAMETERS                            
   !****************************************************************************
   ! Broadcast to all processors. Do not broadcast <verbose>.
   !****************************************************************************
   subroutine broadcast_parameters
      use mp, only: broadcast

      implicit none
      call broadcast(xyzs_layout)
      call broadcast(vms_layout)
      call broadcast(kymus_layout)
      call broadcast(mat_gen)
      call broadcast(mat_read)
      call broadcast(lu_option_switch)
      call broadcast(fields_kxkyz)

   end subroutine broadcast_parameters
    
   !****************************************************************************
   !                     INIT DISTRIBUTION FUNCTION LAYOUTS                     
   !****************************************************************************
   ! Initialize the four layouts that are needed for the distribution functions.
   ! The layouts depend on the grid divisions along (x,y,kx,ky,z,mu,v,species,tubes).
   !**************************************************************************** 
   subroutine init_dist_fn_layouts(nzgrid, ntubes, naky, nakx, nvgrid, nmu, nspec, ny, nx, nalpha)

      use parameters_physics, only: full_flux_surface
      use parameters_physics, only: include_parallel_nonlinearity
      use parameters_numerical, only: split_parallel_dynamics

      implicit none

      integer, intent(in) :: nzgrid, ntubes, naky, nakx, nvgrid, nmu, nspec, ny, nx, nalpha
      integer :: nz, nvpa

      !-------------------------------------------------------------------------
      
      ! Derived grid variables 
      nz = 2*nzgrid + 1
      nvpa = 2*nvgrid
      
      ! Print some basic info to the output file
      call print_info_to_command_prompt
      
      ! Perform sanity checks for the integers
      call sanity_checks

      ! Initialise the parallelisation layouts
      call init_vmu_layout(nzgrid, ntubes, naky, nakx, nvgrid, nmu, nspec, ny, nx, nalpha)
      call init_kxkyz_layout(nzgrid, ntubes, naky, nakx, nvgrid, nmu, nspec) 
      if (full_flux_surface) call init_kxyz_layout(nzgrid, ntubes, naky, nakx, nvgrid, nmu, nspec, ny)
      if (include_parallel_nonlinearity) call init_xyz_layout(nzgrid, ntubes, naky, nakx, nvgrid, nmu, nspec, ny, nx)
      if (.not. split_parallel_dynamics) call init_kymus_layout(nzgrid, ntubes, naky, nakx, nvgrid, nmu, nspec)
   
   contains
   
      !=========================================================================
      ! Sanity checks: on a 32-bit system (4 bytes) the largest signed integer 
      ! is 2,147,483,647. Integers which are bigger will run into overflow errors.
      ! In this case, we would need to define these integers using 64 bits instead.
      !=========================================================================
      subroutine sanity_checks
      
         use mp, only: mp_abort
         use iso_fortran_env, only: int64, int32

         implicit none
         
         integer(int64) :: prod_ivmu, prod_ikxkyz
         integer(int64), parameter :: int32_warn = 1500000000_int64
         integer(int64), parameter :: int32_max = huge(0_int32)

         !----------------------------------------------------------------------
         
         ! Check the ivmu and ikxkyz indices using 64-bit integers
         prod_ivmu = int(nmu, int64) * int(nvpa, int64) * int(nspec, int64)
         prod_ikxkyz = int(nakx, int64) * int(naky, int64) * int(nz, int64) * int(ntubes, int64) * int(nspec, int64)
      
         ! Write warnings if we're getting close to the integer threshold
         if (prod_ivmu > int32_warn) write (*,*) 'WARNING: The index "ivmu" is getting close to 32-bit integer overflow.'
         if (prod_ikxkyz > int32_warn) write (*,*) 'WARNING: The index "ikxkyz" is getting close to 32-bit integer overflow.'

         ! Hard failure if truly exceeding 32-bit range
         if (prod_ivmu > int32_max) call mp_abort('Index "ivmu" exceeds 32-bit integer range; use 64-bit integers instead. Aborting.')
         if (prod_ikxkyz > int32_max) call mp_abort('Index "ikxkyz" exceeds 32-bit integer range; use 64-bit integers. Aborting.')

      end subroutine sanity_checks
      
      !=========================================================================
      ! Print some basic info to the output file 
      ! Note that Gigabyte = 10**9 and Gigibyte = 2**20 = 1024**3
      !=========================================================================
      subroutine print_info_to_command_prompt

         use mp, only: nproc
         use iso_fortran_env, only: int64
         use debug_flags, only: print_extra_info_to_terminal

         implicit none
         
         integer(int64) :: prod_points

         !----------------------------------------------------------------------
         
         ! Only print the number of processors
         if (print_extra_info_to_terminal) then
            write (*, '(A)') '############################################################'
            write (*, '(A)') '                       PARALLELISATION                      '
            write (*, '(A)') '############################################################'
            if (nproc == 1) then
               write (*, *) ''; write (*, '(A,I0,A)') ' Running on ', nproc, ' processor.'
            else
               write (*, *) ''; write (*, '(A,I0,A)') ' Running on ', nproc, ' processors.'
            end if
         end if
         
         ! Also print the grids and the number of points for 
         ! the distribution function and the electrostatic potential
         if (verbose .and. print_extra_info_to_terminal) then
         
            ! Grids
            write (*, *) ''
            write (*, *) 'Grids:'
            write (*, '(A,I0)') '    nx = ', nx
            write (*, '(A,I0)') '    ny = ', ny
            write (*, '(A,I0)') '    nz = ', nz
            write (*, '(A,I0)') '    nkx = ', nakx
            write (*, '(A,I0)') '    nky = ', naky
            write (*, '(A,I0)') '    nmu = ', nmu
            write (*, '(A,I0)') '    nvpa = ', nvpa
            write (*, '(A,I0)') '    nspec = ', nspec
            write (*, '(A,I0)') '    ntubes = ', ntubes
            write (*, '(A,I0)') '    nalpha = ', nalpha
            write (*, *) ''
            
            ! Total number of points (which can be more than the max 32-bit integer)
            write (*, *) 'Total number of grid points:'
            prod_points = int(nx, int64) * int(ny, int64) * int(nz, int64) * int(nmu, int64) * int(nvpa, int64) * int(nspec, int64)
            write (*, '(A,I0,A,F0.2,A,F0.2,A)') '    nx*ny*nz*nmu*nvpa*nspec = ', prod_points, &
               ' = ', prod_points/1024./1024., ' MiB = ', prod_points/1024./1024./1024., ' GiB.'
            prod_points = int(nakx, int64) * int(naky, int64) * int(nz, int64) * int(nmu, int64) * int(nvpa, int64) * int(nspec, int64)
            write (*, '(A,I0,A,F0.2,A,F0.2,A)') '    nkx*nky*nz*nmu*nvpa*nspec = ', prod_points, &
               ' = ', prod_points/1024./1024., ' MiB = ', prod_points/1024./1024./1024., ' GiB.'
            prod_points = int(nakx, int64) * int(naky, int64) * int(nz, int64) * int(ntubes, int64) * int(nspec, int64)
            write (*, '(A,I0,A,F0.2,A,F0.2,A)') '    nkx*nky*nz = ', prod_points, &
               ' = ', prod_points/1024./1024., ' MiB = ', prod_points/1024./1024./1024., ' GiB.'
            write (*, *) ''
            write (*, *) 'Number of points to be parallelised:'
            write (*, '(A,I10,A)') '    vmu-layout:  ', nmu*nvpa*nspec, '    (nmu*nvpa*nspec)'
            write (*, '(A,I10,A)') '    kxkyz-layout:', nakx*naky*nz*ntubes*nspec, '    (nkx*nky*nz*ntubes*nspec)'
            if (full_flux_surface) write (*, '(A,I7,A)') '    kxyz-layout: ', (nakx/2+1)*ny*nz*ntubes*nspec, '    ((nakx/2+1)*ny*nz*ntubes*nspec)'
            if (include_parallel_nonlinearity) write (*, '(A,I7,A)') '    xyz-layout:  ', nx*ny*nz*nspec, '    (nx*ny*nz*nspec)'
            if (.not. split_parallel_dynamics) write (*, '(A,I7,A)') '    kymus-layout:', naky*nmu*nspec, '    (nky*nmu*nspec)'
            write (*, *) ''
            write (*, *) 'Number of points per processor:'
            write (*, '(A,I10)') '    vmu-layout:  ', (nmu*nvpa*nspec-1)/nproc+1
            write (*, '(A,I10)') '    kxkyz-layout:', (nakx*naky*nz*ntubes*nspec-1)/nproc+1
            if (full_flux_surface) write (*, '(A,I7)') '    kxyz-layout: ', ((nakx/2+1)*ny*nz*ntubes*nspec-1)/nproc+1
            if (include_parallel_nonlinearity) write (*, '(A,I7)') '    xyz-layout:  ', (nx*ny*nz*nspec-1)/nproc+1
            if (.not. split_parallel_dynamics) write (*, '(A,I7)') '    kymus-layout:', (naky*nmu*nspec-1)/nproc+1
         end if
      
      end subroutine print_info_to_command_prompt
      
   end subroutine init_dist_fn_layouts
   
!###############################################################################
!                                  KXKYZ LAYOUT                                 
!###############################################################################
! Define the layout/division of the (kx-ky-z) grid on the processors.
! 
! Divide (naky*nakx*nzed*ntubes*nspec - 1) points over the processors.
! Each processor will be in charge of <blocksize> points ranging from <llim_proc>
! up to <ulim_proc>. The "is_idx", "ikx_idx", "iky_idx", "iz_idx" and "it_idx"
! routines will return the values of (species, kx, ky, z, tube) associated with
! point <i> = <idx> = <ikxkyz>.
! 
! The KXKYZ Layout options are "xyzs", "xzys", "yxzs", "zxys", "yzxs", which
! decide the order of (kx, ky, z, species, tubes) of the points <idx>.
!       x --> [1, ..., nakx]
!       y --> [1, ..., naky]
!       z --> [1, ..., nzed = 2*nzgrid+1]
!       s --> [1, ..., nspec]
!       t --> [1, ..., ntubes]
! 
! The "proc_id" routine returns the id of the processor.
!       proc_id --> [0, ..., nproc-1]
! 
! The "idx_kxkyz" routine returns the id of the point.
!       idx --> [0, ..., ulim_world] --> point (iky,ikx,iz,it,is)
! 
! The "idx_local" routine returns true if the point <idx> lies on the current
! processor <iproc> and returns false if the point <idx> does not lie on it.
!###############################################################################

   !****************************************************************************
   !                          INITIALISE KXKYZ LAYOUT                           
   !****************************************************************************
   subroutine init_kxkyz_layout(nzgrid, ntubes, naky, nakx, nvgrid, nmu, nspec)

      use mp, only: iproc, nproc

      implicit none

      integer, intent(in) :: nzgrid, ntubes, naky, nakx, nvgrid, nmu, nspec
      logical, save :: kxkyz_initialised = .false.
      
      !-------------------------------------------------------------------------

      ! Only initialise once
      if (kxkyz_initialised) return
      kxkyz_initialised = .true.

      ! Tubes and species
      kxkyz_lo%nspec = nspec
      kxkyz_lo%ntubes = ntubes
      
      ! Velocity grid
      kxkyz_lo%nmu = nmu
      kxkyz_lo%nvgrid = nvgrid
      kxkyz_lo%nvpa = 2 * nvgrid
      
      ! Parallel space grid
      kxkyz_lo%nzgrid = nzgrid
      kxkyz_lo%nzed = 2 * nzgrid + 1
      
      ! Perpendicular space grid
      kxkyz_lo%naky = naky
      kxkyz_lo%nakx = nakx

      ! The current processor
      kxkyz_lo%iproc = iproc
      
      ! Number of points that need to be divided over the processors
      kxkyz_lo%llim_world = 0
      kxkyz_lo%ulim_world = naky * nakx * kxkyz_lo%nzed * ntubes * nspec - 1
      
      ! Number of points per processor
      kxkyz_lo%blocksize = kxkyz_lo%ulim_world / nproc + 1
      kxkyz_lo%llim_proc = kxkyz_lo%blocksize * iproc
      kxkyz_lo%ulim_proc = min(kxkyz_lo%ulim_world, kxkyz_lo%llim_proc + kxkyz_lo%blocksize - 1)
      kxkyz_lo%ulim_alloc = max(kxkyz_lo%llim_proc, kxkyz_lo%ulim_proc)

   end subroutine init_kxkyz_layout

   !****************************************************************************
   !                                   IS_IDX                                   
   !****************************************************************************
   ! Returns the species index, in the range [1, ..., nspec], of point i.
   ! The order of the division does not matter, so no need for branching
   !****************************************************************************
   elemental function is_idx_kxkyz(lo, i)
   
      implicit none
      
      integer :: is_idx_kxkyz
      type(kxkyz_layout_type), intent(in) :: lo
      integer, intent(in) :: i
      
      is_idx_kxkyz = 1 + mod((i - lo%llim_world) / lo%nakx / lo%naky / lo%nzed / lo%ntubes, lo%nspec)
      
   end function is_idx_kxkyz

   !****************************************************************************
   !                                  IKX_IDX                                   
   !****************************************************************************
   ! Returns the kx index, in the range [1, ..., nakx], of point i.
   !****************************************************************************
   elemental function ikx_idx_kxkyz(lo, i)

      implicit none

      integer :: ikx_idx_kxkyz
      type(kxkyz_layout_type), intent(in) :: lo
      integer, intent(in) :: i

      select case (xyzs_layout)
      case ('xyzs')
         ikx_idx_kxkyz = 1 + mod((i - lo%llim_world), lo%nakx)
      case ('xzys')
         ikx_idx_kxkyz = 1 + mod((i - lo%llim_world), lo%nakx)
      case ('yxzs')
         ikx_idx_kxkyz = 1 + mod((i - lo%llim_world) / lo%naky, lo%nakx)
      case ('zxys')
         ikx_idx_kxkyz = 1 + mod((i - lo%llim_world) / lo%nzed / lo%ntubes, lo%nakx)
      case ('zyxs')
         ikx_idx_kxkyz = 1 + mod((i - lo%llim_world) / lo%nzed / lo%ntubes / lo%naky, lo%nakx)
      case ('yzxs')
         ikx_idx_kxkyz = 1 + mod((i - lo%llim_world) / lo%naky / lo%nzed / lo%ntubes, lo%nakx)
      end select

   end function ikx_idx_kxkyz

   !****************************************************************************
   !                                  IKY_IDX                                   
   !****************************************************************************
   ! Returns the ky index, in the range [1, ..., naky], of point i.
   !****************************************************************************
   elemental function iky_idx_kxkyz(lo, i)

      implicit none
      
      integer :: iky_idx_kxkyz
      type(kxkyz_layout_type), intent(in) :: lo
      integer, intent(in) :: i

      select case (xyzs_layout)
      case ('yxzs')
         iky_idx_kxkyz = 1 + mod(i - lo%llim_world, lo%naky)
      case ('yzxs')
         iky_idx_kxkyz = 1 + mod(i - lo%llim_world, lo%naky)
      case ('xyzs')
         iky_idx_kxkyz = 1 + mod((i - lo%llim_world) / lo%nakx, lo%naky)
      case ('zyxs')
         iky_idx_kxkyz = 1 + mod((i - lo%llim_world) / lo%nzed / lo%ntubes, lo%naky)
      case ('zxys')
         iky_idx_kxkyz = 1 + mod((i - lo%llim_world) / lo%nzed / lo%ntubes / lo%nakx, lo%naky)
      case ('xzys')
         iky_idx_kxkyz = 1 + mod((i - lo%llim_world) / lo%nakx / lo%nzed / lo%ntubes, lo%naky)
      end select

   end function iky_idx_kxkyz

   !****************************************************************************
   !                                   IZ_IDX                                   
   !****************************************************************************
   ! Returns the z index, in the range [1, ..., nzed], of point i.
   !****************************************************************************
   elemental function iz_idx_kxkyz(lo, i)

      implicit none
      
      integer :: iz_idx_kxkyz
      type(kxkyz_layout_type), intent(in) :: lo
      integer, intent(in) :: i

      select case (xyzs_layout)
      case ('zyxs')
         iz_idx_kxkyz = -lo%nzgrid + mod((i - lo%llim_world), lo%nzed)
      case ('zxys')
         iz_idx_kxkyz = -lo%nzgrid + mod((i - lo%llim_world), lo%nzed)
      case ('xzys')
         iz_idx_kxkyz = -lo%nzgrid + mod((i - lo%llim_world) / lo%nakx, lo%nzed)
      case ('yzxs')
         iz_idx_kxkyz = -lo%nzgrid + mod((i - lo%llim_world) / lo%naky, lo%nzed)
      case ('yxzs')
         iz_idx_kxkyz = -lo%nzgrid + mod((i - lo%llim_world) / lo%naky / lo%nakx, lo%nzed)
      case ('xyzs')
         iz_idx_kxkyz = -lo%nzgrid + mod((i - lo%llim_world) / lo%nakx / lo%naky, lo%nzed)
      end select

   end function iz_idx_kxkyz
   
   !****************************************************************************
   !                                   IT_IDX                                   
   !****************************************************************************
   ! Returns the tubes index, in the range [1, ..., ntubes], of point i.
   !****************************************************************************
   elemental function it_idx_kxkyz(lo, i)

      implicit none
      
      integer :: it_idx_kxkyz
      type(kxkyz_layout_type), intent(in) :: lo
      integer, intent(in) :: i

      select case (xyzs_layout)
      case ('zyxs')
         it_idx_kxkyz = 1 + mod((i - lo%llim_world) / lo%nzed, lo%ntubes)
      case ('zxys')
         it_idx_kxkyz = 1 + mod((i - lo%llim_world) / lo%nzed, lo%ntubes)
      case ('xzys')
         it_idx_kxkyz = 1 + mod((i - lo%llim_world) / lo%nzed / lo%nakx, lo%ntubes)
      case ('yzxs')
         it_idx_kxkyz = 1 + mod((i - lo%llim_world) / lo%nzed / lo%naky, lo%ntubes)
      case ('yxzs')
         it_idx_kxkyz = 1 + mod((i - lo%llim_world) / lo%nzed / lo%naky / lo%nakx, lo%ntubes)
      case ('xyzs')
         it_idx_kxkyz = 1 + mod((i - lo%llim_world) / lo%nzed / lo%nakx / lo%naky, lo%ntubes)
      end select

   end function it_idx_kxkyz

   !****************************************************************************
   !                                   PROC_ID                                  
   !****************************************************************************
   ! Returns the processor index, in the range [0, ..., nproc-1], of point i.
   !****************************************************************************
   elemental function proc_id_kxkyz(lo, i)
   
      implicit none
      
      integer :: proc_id_kxkyz
      type(kxkyz_layout_type), intent(in) :: lo
      integer, intent(in) :: i

      proc_id_kxkyz = i / lo%blocksize

   end function proc_id_kxkyz

   !****************************************************************************
   !                                  IDX_KXKYZ                                 
   !****************************************************************************
   ! Returns the index <idx> in [0, ..., ulim_world] of point (iky,ikx,iz,it,is).
   !****************************************************************************
   elemental function idx_kxkyz(lo, iky, ikx, iz, it, is)

      implicit none

      integer :: idx_kxkyz
      type(kxkyz_layout_type), intent(in) :: lo
      integer, intent(in) :: iky, ikx, iz, it, is

      select case (xyzs_layout)
      case ('xyzs')
         idx_kxkyz = ikx - 1 + lo%nakx * (iky - 1 + lo%naky * (iz + lo%nzgrid + lo%nzed * (it - 1 + lo%ntubes * (is - 1))))
      case ('xzys')
         idx_kxkyz = ikx - 1 + lo%nakx * (iz + lo%nzgrid + lo%nzed * (it - 1 + lo%ntubes * (iky - 1 + lo%naky * (is - 1))))
      case ('yxzs')
         idx_kxkyz = iky - 1 + lo%naky * (ikx - 1 + lo%nakx * (iz + lo%nzgrid + lo%nzed * (it - 1 + lo%ntubes * (is - 1))))
      case ('yzxs')
         idx_kxkyz = iky - 1 + lo%naky * (iz + lo%nzgrid + lo%nzed * (it - 1 + lo%ntubes * (ikx - 1 + lo%nakx * (is - 1))))
      case ('zyxs')
         idx_kxkyz = iz + lo%nzgrid + lo%nzed * (it - 1 + lo%ntubes * (iky - 1 + lo%naky * (ikx - 1 + lo%nakx * (is - 1))))
      case ('zxys')
         idx_kxkyz = iz + lo%nzgrid + lo%nzed * (it - 1 + lo%ntubes * (ikx - 1 + lo%nakx * (iky - 1 + lo%naky * (is - 1))))
      end select

   end function idx_kxkyz

   !****************************************************************************
   !                                 IDX_LOCAL                                  
   !****************************************************************************
   ! When <lo> = <kxkyz_layout_type> then the interface "idx_local" will call
   ! idx_local_kxkyz, next "idx_local" will call iz_local_kxkyz(lo, idx).
   ! Returns true if the point <idx> lies on the current processor <iproc>.
   ! Returns false if the point <idx> does not lie on the current processor <iproc>.
   !****************************************************************************
   elemental function idx_local_kxkyz(lo, iky, ikx, iz, it, is)

      implicit none
      
      logical :: idx_local_kxkyz
      type(kxkyz_layout_type), intent(in) :: lo
      integer, intent(in) :: iky, ikx, iz, it, is

      idx_local_kxkyz = idx_local(lo, idx(lo, iky, ikx, iz, it, is))
      
   end function idx_local_kxkyz

   !****************************************************************************
   !                                 IZ_LOCAL                                   
   !****************************************************************************
   ! Returns true if the point <idx=iz> lies on the current processor <iproc>.
   ! Returns false if the point <idx=iz> does not lie on the current processor <iproc>.
   !****************************************************************************
   elemental function iz_local_kxkyz(lo, iz)
   
      implicit none
      
      logical :: iz_local_kxkyz
      type(kxkyz_layout_type), intent(in) :: lo
      integer, intent(in) :: iz

      iz_local_kxkyz = lo%iproc == proc_id(lo, iz)
      
   end function iz_local_kxkyz

!###############################################################################
!                                  KXYZ LAYOUT                                  
!###############################################################################
! Define the layout/division of the (kx-y-z) grid on the processors.
! 
! Divide (ny*nakx*nzed*ntubes*nspec - 1) points over the processors.
! Each processor will be in charge of <blocksize> points ranging from <llim_proc>
! up to <ulim_proc>. The "is_idx", "ikx_idx", "iy_idx", "iz_idx" and "it_idx"
! routines will return the values of (species, kx, y, z, tube) associated with
! point <i> = <idx> = <ikxyz>.
!
! The KXYZ Layout options are "xyzs", "xzys", "yxzs", "zxys", "yzxs", which
! decide the order of (kx, y, z, species, tubes) of the points <idx>.
!       x --> [1, ..., nakx]
!       y --> [1, ..., ny]
!       z --> [1, ..., nzed = 2*nzgrid+1]
!       s --> [1, ..., nspec]
!       t --> [1, ..., ntubes]
! 
! The "proc_id" routine returns the id of the processor.
!       proc_id --> [0, ..., nproc-1]
! 
! The "idx_kxkyz" routine returns the id of the point.
!       idx --> [0, ..., ulim_world] --> point (iky,ikx,iz,it,is)
! 
! The "idx_local" routine returns true if the point <idx> lies on the current
! processor <iproc> and returns false if the point <idx> does not lie on it.
!###############################################################################

   !****************************************************************************
   !                         INITITIALIZE KXYZ LAYOUT                           
   !****************************************************************************
   subroutine init_kxyz_layout(nzgrid, ntubes, naky, nakx, nvgrid, nmu, nspec, ny)

      use mp, only: iproc, nproc

      implicit none

      integer, intent(in) :: nzgrid, ntubes, ny, naky, nakx, nvgrid, nmu, nspec
      logical, save :: kxyz_initialised = .false.
      
      !-------------------------------------------------------------------------

      ! Only initialise once
      if (kxyz_initialised) return
      kxyz_initialised = .true.

      ! Tubes and species
      kxyz_lo%nspec = nspec
      kxyz_lo%ntubes = ntubes

      ! Velocity grid
      kxyz_lo%nmu = nmu
      kxyz_lo%nvgrid = nvgrid
      kxyz_lo%nvpa = 2*nvgrid

      ! Parallel space grid
      kxyz_lo%nzgrid = nzgrid
      kxyz_lo%nzed = 2*nzgrid+1

      ! Perpendicular space grid
      kxyz_lo%ny = ny
      kxyz_lo%naky = naky
      kxyz_lo%nakx = nakx
      kxyz_lo%ikx_max = nakx/2 + 1

      ! The current processor
      kxyz_lo%iproc = iproc

      ! Number of points that need to be divided over the processors
      kxyz_lo%llim_world = 0
      kxyz_lo%ulim_world = ny * kxyz_lo%ikx_max * kxyz_lo%nzed * ntubes * nspec - 1

      ! Number of points per processor
      kxyz_lo%blocksize = kxyz_lo%ulim_world/nproc + 1
      kxyz_lo%llim_proc = kxyz_lo%blocksize*iproc
      kxyz_lo%ulim_proc = min(kxyz_lo%ulim_world, kxyz_lo%llim_proc + kxyz_lo%blocksize - 1)
      kxyz_lo%ulim_alloc = max(kxyz_lo%llim_proc, kxyz_lo%ulim_proc)

   end subroutine init_kxyz_layout

   !****************************************************************************
   !                                   IS_IDX                                   
   !****************************************************************************
   ! Returns the species index, in the range [1, ..., nspec], of point i.
   ! The order of the division does not matter, so no need for branching
   !****************************************************************************
   elemental function is_idx_kxyz(lo, i)
   
      implicit none
      
      integer :: is_idx_kxyz
      type(kxyz_layout_type), intent(in) :: lo
      integer, intent(in) :: i
      
      is_idx_kxyz = 1 + mod((i - lo%llim_world) / lo%ikx_max / lo%ny / lo%nzed / lo%ntubes, lo%nspec)
      
   end function is_idx_kxyz

   !****************************************************************************
   !                                  IKX_IDX                                   
   !****************************************************************************
   ! Returns the kx index, in the range [1, ..., nakx], of point i.
   !****************************************************************************
   elemental function ikx_idx_kxyz(lo, i)

      implicit none

      integer :: ikx_idx_kxyz
      type(kxyz_layout_type), intent(in) :: lo
      integer, intent(in) :: i

      select case (xyzs_layout)
      case ('xyzs')
         ikx_idx_kxyz = 1 + mod((i - lo%llim_world), lo%ikx_max)
      case ('xzys')
         ikx_idx_kxyz = 1 + mod((i - lo%llim_world), lo%ikx_max)
      case ('yxzs')
         ikx_idx_kxyz = 1 + mod((i - lo%llim_world) / lo%ny, lo%ikx_max)
      case ('zxys')
         ikx_idx_kxyz = 1 + mod((i - lo%llim_world) / lo%nzed / lo%ntubes, lo%ikx_max)
      case ('zyxs')
         ikx_idx_kxyz = 1 + mod((i - lo%llim_world) / lo%nzed / lo%ntubes / lo%ny, lo%ikx_max)
      case ('yzxs')
         ikx_idx_kxyz = 1 + mod((i - lo%llim_world) / lo%ny / lo%nzed / lo%ntubes, lo%ikx_max)
      end select

   end function ikx_idx_kxyz

   !****************************************************************************
   !                                   IY_IDX                                   
   !****************************************************************************
   ! Returns the y index, in the range [1, ..., ny], of point i.
   !****************************************************************************
   elemental function iy_idx_kxyz(lo, i)

      implicit none
      
      integer :: iy_idx_kxyz
      type(kxyz_layout_type), intent(in) :: lo
      integer, intent(in) :: i

      select case (xyzs_layout)
      case ('yxzs')
         iy_idx_kxyz = 1 + mod(i - lo%llim_world, lo%ny)
      case ('yzxs')
         iy_idx_kxyz = 1 + mod(i - lo%llim_world, lo%ny)
      case ('xyzs')
         iy_idx_kxyz = 1 + mod((i - lo%llim_world) / lo%ikx_max, lo%ny)
      case ('zyxs')
         iy_idx_kxyz = 1 + mod((i - lo%llim_world) / lo%nzed / lo%ntubes, lo%ny)
      case ('zxys')
         iy_idx_kxyz = 1 + mod((i - lo%llim_world) / lo%nzed / lo%ntubes / lo%ikx_max, lo%ny)
      case ('xzys')
         iy_idx_kxyz = 1 + mod((i - lo%llim_world) / lo%ikx_max / lo%nzed / lo%ntubes, lo%ny)
      end select

   end function iy_idx_kxyz

   !****************************************************************************
   !                                   IZ_IDX                                   
   !****************************************************************************
   ! Returns the z index, in the range [1, ..., nzed], of point i.
   !****************************************************************************
   elemental function iz_idx_kxyz(lo, i)

      implicit none
      
      integer :: iz_idx_kxyz
      type(kxyz_layout_type), intent(in) :: lo
      integer, intent(in) :: i

      select case (xyzs_layout)
      case ('zyxs')
         iz_idx_kxyz = -lo%nzgrid + mod((i - lo%llim_world), lo%nzed)
      case ('zxys')
         iz_idx_kxyz = -lo%nzgrid + mod((i - lo%llim_world), lo%nzed)
      case ('yzxs')
         iz_idx_kxyz = -lo%nzgrid + mod((i - lo%llim_world) / lo%ny, lo%nzed)
      case ('xzys')
         iz_idx_kxyz = -lo%nzgrid + mod((i - lo%llim_world) / lo%ikx_max, lo%nzed)
      case ('yxzs')
         iz_idx_kxyz = -lo%nzgrid + mod((i - lo%llim_world) / lo%ny / lo%ikx_max, lo%nzed)
      case ('xyzs')
         iz_idx_kxyz = -lo%nzgrid + mod((i - lo%llim_world) / lo%ikx_max / lo%ny, lo%nzed)
      end select

   end function iz_idx_kxyz

   !****************************************************************************
   !                                   IT_IDX                                   
   !****************************************************************************
   ! Returns the tubes index, in the range [1, ..., ntubes], of point i.
   !****************************************************************************
   elemental function it_idx_kxyz(lo, i)

      implicit none
      
      integer :: it_idx_kxyz
      type(kxyz_layout_type), intent(in) :: lo
      integer, intent(in) :: i

      select case (xyzs_layout)
      case ('zyxs')
         it_idx_kxyz = 1 + mod((i - lo%llim_world) / lo%nzed, lo%ntubes)
      case ('zxys')
         it_idx_kxyz = 1 + mod((i - lo%llim_world) / lo%nzed, lo%ntubes)
      case ('yzxs')
         it_idx_kxyz = 1 + mod((i - lo%llim_world) / lo%nzed / lo%ny, lo%ntubes)
      case ('xzys')
         it_idx_kxyz = 1 + mod((i - lo%llim_world) / lo%nzed / lo%ikx_max, lo%ntubes)
      case ('yxzs')
         it_idx_kxyz = 1 + mod((i - lo%llim_world) / lo%nzed / lo%ny / lo%ikx_max, lo%ntubes)
      case ('xyzs')
         it_idx_kxyz = 1 + mod((i - lo%llim_world) / lo%nzed / lo%ikx_max / lo%ny, lo%ntubes)
      end select

   end function it_idx_kxyz

   !****************************************************************************
   !                                   PROC_ID                                  
   !****************************************************************************
   ! Returns the processor index, in the range [0, ..., nproc-1], of point i.
   !****************************************************************************
   elemental function proc_id_kxyz(lo, i)
   
      implicit none
      
      integer :: proc_id_kxyz
      type(kxyz_layout_type), intent(in) :: lo
      integer, intent(in) :: i

      proc_id_kxyz = i / lo%blocksize

   end function proc_id_kxyz

   !****************************************************************************
   !                                  IDX_KXYZ                                 
   !****************************************************************************
   ! Returns the index <idx> in [0, ..., ulim_world] of point (iy,ikx,iz,it,is).
   !****************************************************************************
   elemental function idx_kxyz(lo, iy, ikx, iz, it, is)

      implicit none

      integer :: idx_kxyz
      type(kxyz_layout_type), intent(in) :: lo
      integer, intent(in) :: iy, ikx, iz, it, is

      select case (xyzs_layout)
      case ('xyzs')
         idx_kxyz = ikx - 1 + lo%ikx_max * (iy - 1 + lo%ny * (iz + lo%nzgrid + lo%nzed * (it - 1 + lo%ntubes * (is - 1))))
      case ('xzys')
         idx_kxyz = ikx - 1 + lo%ikx_max * (iz + lo%nzgrid + lo%nzed * (it - 1 + lo%ntubes * (iy - 1 + lo%ny * (is - 1))))
      case ('yxzs')
         idx_kxyz = iy - 1 + lo%ny * (ikx - 1 + lo%ikx_max * (iz + lo%nzgrid + lo%nzed * (it - 1 + lo%ntubes * (is - 1))))
      case ('yzxs')
         idx_kxyz = iy - 1 + lo%ny * (iz + lo%nzgrid + lo%nzed * (it - 1 + lo%ntubes * (ikx - 1 + lo%ikx_max * (is - 1))))
      case ('zyxs')
         idx_kxyz = iz + lo%nzgrid + lo%nzed * (it - 1 + lo%ntubes * (iy - 1 + lo%ny * (ikx - 1)))
      case ('zxys')
         idx_kxyz = iz + lo%nzgrid + lo%nzed * (it - 1 + lo%ntubes * (ikx - 1 + lo%ikx_max * (iy - 1)))
      end select

   end function idx_kxyz

   !****************************************************************************
   !                                 IDX_LOCAL                                  
   !****************************************************************************
   ! When <lo> = <kxyz_layout_type> then the interface "idx_local" will call
   ! idx_local_kxkyz, next "idx_local" will call iz_local_kxyz(lo, idx).
   ! Returns true if the point <idx> lies on the current processor <iproc>.
   ! Returns false if the point <idx> does not lie on the current processor <iproc>.
   !****************************************************************************
   elemental function idx_local_kxyz(lo, iy, ikx, iz, it, is)

      implicit none
      
      logical :: idx_local_kxyz
      type(kxyz_layout_type), intent(in) :: lo
      integer, intent(in) :: iy, ikx, iz, it, is

      idx_local_kxyz = idx_local(lo, idx(lo, iy, ikx, iz, it, is))
      
   end function idx_local_kxyz

   !****************************************************************************
   !                                 IZ_LOCAL                                   
   !****************************************************************************
   ! Returns true if the point <idx=iz> lies on the current processor <iproc>.
   ! Returns false if the point <idx=iz> does not lie on the current processor <iproc>.
   !****************************************************************************
   elemental function iz_local_kxyz(lo, iz)
   
      implicit none
      
      logical :: iz_local_kxyz
      type(kxyz_layout_type), intent(in) :: lo
      integer, intent(in) :: iz

      iz_local_kxyz = lo%iproc == proc_id(lo, iz)
      
   end function iz_local_kxyz

!###############################################################################
!                                   XYZ LAYOUT                                  
!###############################################################################
! Define the layout/division of the (x-y-z) grid on the processors.
! 
! Divide (ny*nx*nzed*ntubes*nspec - 1) points over the processors.
! Each processor will be in charge of <blocksize> points ranging from <llim_proc>
! up to <ulim_proc>. The "is_idx", "ix_idx", "iy_idx", "iz_idx" and "it_idx"
! routines will return the values of (species, kx, y, z, tube) associated with
! point <i> = <idx> = <ikxyz>.
! 
! The KXYZ Layout options are "xyzs", "xzys", "yxzs", "zxys", "yzxs", which
! decide the order of (kx, y, z, species, tubes) of the points <idx>.
!       x --> [1, ..., nx]
!       y --> [1, ..., ny]
!       z --> [1, ..., nzed = 2*nzgrid+1]
!       s --> [1, ..., nspec]
!       t --> [1, ..., ntubes]
! 
! The "proc_id" routine returns the id of the processor.
!       proc_id --> [0, ..., nproc-1]
! 
! The "idx_kxkyz" routine returns the id of the point.
!       idx --> [0, ..., ulim_world] --> point (iky,ikx,iz,it,is)
! 
! The "idx_local" routine returns true if the point <idx> lies on the current
! processor <iproc> and returns false if the point <idx> does not lie on it.
!###############################################################################

   subroutine init_xyz_layout(nzgrid, ntubes, naky, nakx, nvgrid, nmu, nspec, ny, nx)

      use mp, only: iproc, nproc

      implicit none

      integer, intent(in) :: nzgrid, ntubes, ny, nx, naky, nakx, nvgrid, nmu, nspec
      logical, save :: xyz_initialised = .false.
      
      !-------------------------------------------------------------------------

      ! Only initialise once
      if (xyz_initialised) return
      xyz_initialised = .true.
   
      ! Tubes and species
      xyz_lo%nspec = nspec
      xyz_lo%ntubes = ntubes

      ! Velocity grid
      xyz_lo%nmu = nmu
      xyz_lo%nvgrid = nvgrid
      xyz_lo%nvpa = 2*nvgrid

      ! Parallel space grid
      xyz_lo%nzgrid = nzgrid
      xyz_lo%nzed = 2*nzgrid+1

      ! Perpendicular space grid
      xyz_lo%ny = ny
      xyz_lo%nx = nx
      xyz_lo%naky = naky
      xyz_lo%nakx = nakx

      ! The current processor
      xyz_lo%iproc = iproc

      ! Number of points that need to be divided over the processors
      xyz_lo%llim_world = 0
      xyz_lo%ulim_world = ny*nx*xyz_lo%nzed*xyz_lo%ntubes*nspec - 1

      ! Number of points per processor
      xyz_lo%blocksize = xyz_lo%ulim_world/nproc + 1
      xyz_lo%llim_proc = xyz_lo%blocksize*iproc
      xyz_lo%ulim_proc = min(xyz_lo%ulim_world, xyz_lo%llim_proc + xyz_lo%blocksize - 1)
      xyz_lo%ulim_alloc = max(xyz_lo%llim_proc, xyz_lo%ulim_proc)

   end subroutine init_xyz_layout

   !****************************************************************************
   !                                   IS_IDX                                   
   !****************************************************************************
   ! Returns the species index, in the range [1, ..., nspec], of point i.
   ! The order of the division does not matter, so no need for branching
   !****************************************************************************
   elemental function is_idx_xyz(lo, i)
      implicit none
      
      integer :: is_idx_xyz
      type(xyz_layout_type), intent(in) :: lo
      integer, intent(in) :: i
      
      is_idx_xyz = 1 + mod((i - lo%llim_world) / lo%nx / lo%ny / lo%nzed / lo%ntubes, lo%nspec)
      
   end function is_idx_xyz

   !****************************************************************************
   !                                   IX_IDX                                   
   !****************************************************************************
   ! Returns the x index, in the range [1, ..., nx], of point i.
   !****************************************************************************
   elemental function ix_idx_xyz(lo, i)

      implicit none

      integer :: ix_idx_xyz
      type(xyz_layout_type), intent(in) :: lo
      integer, intent(in) :: i

      select case (xyzs_layout)
      case ('xyzs')
         ix_idx_xyz = 1 + mod((i - lo%llim_world), lo%nx)
      case ('xzys')
         ix_idx_xyz = 1 + mod((i - lo%llim_world), lo%nx)
      case ('yxzs')
         ix_idx_xyz = 1 + mod((i - lo%llim_world) / lo%ny, lo%nx)
      case ('zxys')
         ix_idx_xyz = 1 + mod((i - lo%llim_world) / lo%nzed / lo%ntubes, lo%nx)
      case ('zyxs')
         ix_idx_xyz = 1 + mod((i - lo%llim_world) / lo%nzed / lo%ntubes / lo%ny, lo%nx)
      case ('yzxs')
         ix_idx_xyz = 1 + mod((i - lo%llim_world) / lo%ny / lo%nzed / lo%ntubes, lo%nx)
      end select

   end function ix_idx_xyz

   !****************************************************************************
   !                                   IY_IDX                                   
   !****************************************************************************
   ! Returns the y index, in the range [1, ..., ny], of point i.
   !****************************************************************************
   elemental function iy_idx_xyz(lo, i)

      implicit none
      
      integer :: iy_idx_xyz
      type(xyz_layout_type), intent(in) :: lo
      integer, intent(in) :: i

      select case (xyzs_layout)
      case ('yxzs')
         iy_idx_xyz = 1 + mod(i - lo%llim_world, lo%ny)
      case ('yzxs')
         iy_idx_xyz = 1 + mod(i - lo%llim_world, lo%ny)
      case ('xyzs')
         iy_idx_xyz = 1 + mod((i - lo%llim_world) / lo%nx, lo%ny)
      case ('zyxs')
         iy_idx_xyz = 1 + mod((i - lo%llim_world) / lo%nzed / lo%ntubes, lo%ny)
      case ('zxys')
         iy_idx_xyz = 1 + mod((i - lo%llim_world) / lo%nzed / lo%ntubes / lo%nx, lo%ny)
      case ('xzys')
         iy_idx_xyz = 1 + mod((i - lo%llim_world) / lo%nx / lo%nzed / lo%ntubes, lo%ny)
      end select

   end function iy_idx_xyz

   !****************************************************************************
   !                                   IZ_IDX                                   
   !****************************************************************************
   ! Returns the z index, in the range [1, ..., nzed], of point i.
   !****************************************************************************
   elemental function iz_idx_xyz(lo, i)

      implicit none
      
      integer :: iz_idx_xyz
      type(xyz_layout_type), intent(in) :: lo
      integer, intent(in) :: i

      select case (xyzs_layout)
      case ('zyxs')
         iz_idx_xyz = -lo%nzgrid + mod((i - lo%llim_world), lo%nzed)
      case ('zxys')
         iz_idx_xyz = -lo%nzgrid + mod((i - lo%llim_world), lo%nzed)
      case ('yzxs')
         iz_idx_xyz = -lo%nzgrid + mod((i - lo%llim_world) / lo%ny, lo%nzed)
      case ('xzys')
         iz_idx_xyz = -lo%nzgrid + mod((i - lo%llim_world) / lo%nx, lo%nzed)
      case ('yxzs')
         iz_idx_xyz = -lo%nzgrid + mod((i - lo%llim_world) / lo%ny / lo%nx, lo%nzed)
      case ('xyzs')
         iz_idx_xyz = -lo%nzgrid + mod((i - lo%llim_world) / lo%nx / lo%ny, lo%nzed)
      end select

   end function iz_idx_xyz

   !****************************************************************************
   !                                   IT_IDX                                   
   !****************************************************************************
   ! Returns the tubes index, in the range [1, ..., ntubes], of point i.
   !****************************************************************************
   elemental function it_idx_xyz(lo, i)

      implicit none
      
      integer :: it_idx_xyz
      type(xyz_layout_type), intent(in) :: lo
      integer, intent(in) :: i

      select case (xyzs_layout)
      case ('zyxs')
         it_idx_xyz = 1 + mod((i - lo%llim_world) / lo%nzed, lo%ntubes)
      case ('zxys')
         it_idx_xyz = 1 + mod((i - lo%llim_world) / lo%nzed, lo%ntubes)
      case ('yzxs')
         it_idx_xyz = 1 + mod((i - lo%llim_world) / lo%nzed / lo%ny, lo%ntubes)
      case ('xzys')
         it_idx_xyz = 1 + mod((i - lo%llim_world) / lo%nzed / lo%nx, lo%ntubes)
      case ('yxzs')
         it_idx_xyz = 1 + mod((i - lo%llim_world) / lo%nzed / lo%ny / lo%nx, lo%ntubes)
      case ('xyzs')
         it_idx_xyz = 1 + mod((i - lo%llim_world) / lo%nzed / lo%nx / lo%ny, lo%ntubes)
      end select

   end function it_idx_xyz

   !****************************************************************************
   !                                   PROC_ID                                  
   !****************************************************************************
   ! Returns the processor index, in the range [0, ..., nproc-1], of point i.
   !****************************************************************************
   elemental function proc_id_xyz(lo, i)
   
      implicit none
      
      integer :: proc_id_xyz
      type(xyz_layout_type), intent(in) :: lo
      integer, intent(in) :: i

      proc_id_xyz = i / lo%blocksize

   end function proc_id_xyz

   !****************************************************************************
   !                                   IDX_XYZ                                  
   !****************************************************************************
   ! Returns the index <idx> in [0, ..., ulim_world] of point (iy,ix,iz,it,is).
   !****************************************************************************
   elemental function idx_xyz(lo, iy, ix, iz, it, is)

      implicit none

      integer :: idx_xyz
      type(xyz_layout_type), intent(in) :: lo
      integer, intent(in) :: iy, ix, iz, it, is

      select case (xyzs_layout)
      case ('xyzs')
         idx_xyz = ix - 1 + lo%nx * (iy - 1 + lo%ny * (iz + lo%nzgrid + lo%nzed * (it - 1 + lo%ntubes * (is - 1))))
      case ('xzys')
         idx_xyz = ix - 1 + lo%nx * (iz + lo%nzgrid + lo%nzed * (it - 1 + lo%ntubes * (iy - 1 + lo%ny * (is - 1))))
      case ('yxzs')
         idx_xyz = iy - 1 + lo%ny * (ix - 1 + lo%nx * (iz + lo%nzgrid + lo%nzed * (it - 1 + lo%ntubes * (is - 1))))
      case ('yzxs')
         idx_xyz = iy - 1 + lo%ny * (iz + lo%nzgrid + lo%nzed * (it - 1 + lo%ntubes * (ix - 1 + lo%nx * (is - 1))))
      case ('zyxs')
         idx_xyz = iz + lo%nzgrid + lo%nzed * (it - 1 + lo%ntubes * (iy - 1 + lo%ny * (ix - 1)))
      case ('zxys')
         idx_xyz = iz + lo%nzgrid + lo%nzed * (it - 1 + lo%ntubes * (ix - 1 + lo%nx * (iy - 1)))
      end select

   end function idx_xyz

   !****************************************************************************
   !                                 IDX_LOCAL                                  
   !****************************************************************************
   ! When <lo> = <xyz_layout_type> then the interface "idx_local" will call
   ! idx_local_kxkyz, next "idx_local" will call iz_local_kxyz(lo, idx).
   ! Returns true if the point <idx> lies on the current processor <iproc>.
   ! Returns false if the point <idx> does not lie on the current processor <iproc>.
   !****************************************************************************
   elemental function idx_local_xyz(lo, iy, ix, iz, it, is)

      implicit none
      
      logical :: idx_local_xyz
      type(xyz_layout_type), intent(in) :: lo
      integer, intent(in) :: iy, ix, iz, it, is

      idx_local_xyz = idx_local(lo, idx(lo, iy, ix, iz, it, is))
      
   end function idx_local_xyz

   !****************************************************************************
   !                                  IZ_LOCAL                                  
   !****************************************************************************
   ! Returns true if the point <idx=iz> lies on the current processor <iproc>.
   ! Returns false if the point <idx=iz> does not lie on the current processor <iproc>.
   !****************************************************************************
   elemental function iz_local_xyz(lo, iz)
   
      implicit none
      
      logical :: iz_local_xyz
      type(xyz_layout_type), intent(in) :: lo
      integer, intent(in) :: iz

      iz_local_xyz = lo%iproc == proc_id(lo, iz)
      
   end function iz_local_xyz

!###############################################################################
!                                  VMU LAYOUT                                   
!###############################################################################
! Define the layout/division of the (v, mu, s) grid on the processors.
! 
! Divide (2*nvgrid*nmu*nspec - 1) points over the total amount of processors.
! Each processor will be in charge of <blocksize> points ranging from <llim_proc>
! up to <ulim_proc>. The "is_idx", "imu_idx", "iv_idx" routines will return the
! values of (species, mu, vgrid) associated with point <i>=<idx>=<ivmu>.
! 
! The VMU Layout options are "vms" and "mvs" which decide the order of the points.
! First iterate over the first letter, with constant second and third letter, ...
!       v --> [1, ..., 2*nvgrid = nvpa]
!       m --> [1, ..., nmu]
!       s --> [1, ..., nspec]
! 
! The "proc_id" routine returns the id of the processor.
!       proc_id --> [0, ..., nproc-1]
! 
! The "idx_vmu" routine returns the id of the point.
!       idx --> [0, ..., ulim_world] --> point(iv,imu,is)
! 
! The "idx_local" routine returns true if the point <idx> lies on the current
! processor <iproc> and returns false if the point <idx> does not lie on it.
!###############################################################################

   !****************************************************************************
   !                          INITITIALIZE VMU LAYOUT                           
   !****************************************************************************
   subroutine init_vmu_layout(nzgrid, ntubes, naky, nakx, nvgrid, nmu, nspec, ny, nx, nalpha)

      use mp, only: iproc, nproc

      implicit none

      integer, intent(in) :: nzgrid, ntubes, naky, nakx, nvgrid, nmu, nspec, ny, nx, nalpha
      logical, save :: vmu_initialised = .false.
      
      !-------------------------------------------------------------------------

      ! Only initialise once
      if (vmu_initialised) return
      vmu_initialised = .true.
   
      ! Layout, tubes and species
      vmu_lo%xyz = .true.
      vmu_lo%ntubes = ntubes
      vmu_lo%nspec = nspec

      ! Velocity grid
      vmu_lo%nmu = nmu
      vmu_lo%nvgrid = nvgrid
      vmu_lo%nvpa = 2*nvgrid

      ! Parallel space grid
      vmu_lo%nzgrid = nzgrid
      vmu_lo%nzed = 2*nzgrid+1

      ! Perpendicular space grid
      vmu_lo%nx = nx
      vmu_lo%ny = ny
      vmu_lo%nakx = nakx
      vmu_lo%naky = naky
      vmu_lo%nalpha = nalpha

      ! The current processor
      vmu_lo%iproc = iproc

      ! Number of points that need to be divided over the processors
      vmu_lo%llim_world = 0
      vmu_lo%ulim_world = vmu_lo%nvpa*nmu*nspec - 1

      ! Number of points per processor
      vmu_lo%blocksize = vmu_lo%ulim_world/nproc + 1
      vmu_lo%llim_proc = vmu_lo%blocksize*iproc
      vmu_lo%ulim_proc = min(vmu_lo%ulim_world, vmu_lo%llim_proc + vmu_lo%blocksize - 1)
      vmu_lo%ulim_alloc = max(vmu_lo%llim_proc, vmu_lo%ulim_proc)

   end subroutine init_vmu_layout

   !****************************************************************************
   !                                   IS_IDX                                   
   !****************************************************************************
   ! Returns the species index, in the range [1, ..., nspec], of point i.
   ! The order of the division does not matter, so no need for branching
   !****************************************************************************
   elemental function is_idx_vmu(lo, i)

      implicit none
      
      integer :: is_idx_vmu
      type(vmu_layout_type), intent(in) :: lo
      integer, intent(in) :: i

      is_idx_vmu = 1 + mod((i - lo%llim_world) / lo%nvpa / lo%nmu, lo%nspec)

   end function is_idx_vmu

   !****************************************************************************
   !                                  IMU_IDX                                   
   !****************************************************************************
   ! Returns the mu index, in the range [1, ..., nmu], of point i.
   !****************************************************************************
   elemental function imu_idx_vmu(lo, i)

      implicit none

      integer :: imu_idx_vmu
      type(vmu_layout_type), intent(in) :: lo
      integer, intent(in) :: i

      select case (vms_layout)
      case ('vms')
         imu_idx_vmu = 1 + mod((i - lo%llim_world) / lo%nvpa, lo%nmu)
      case ('mvs')
         imu_idx_vmu = 1 + mod((i - lo%llim_world), lo%nmu)
      end select

   end function imu_idx_vmu

   !****************************************************************************
   !                                   IV_IDX                                   
   !****************************************************************************
   ! Returns the v index, in the range [1, ..., 2*nvgrid = nvpa], of point i.
   !****************************************************************************
   elemental function iv_idx_vmu(lo, i)

      implicit none
      
      integer :: iv_idx_vmu
      type(vmu_layout_type), intent(in) :: lo
      integer, intent(in) :: i

      select case (vms_layout)
      case ('vms')
         iv_idx_vmu = 1 + mod((i - lo%llim_world), lo%nvpa)
      case ('mvs')
         iv_idx_vmu = 1 + mod((i - lo%llim_world) / lo%nmu, lo%nvpa)
      end select

   end function iv_idx_vmu

   !****************************************************************************
   !                                   PROC_ID                                  
   !****************************************************************************
   ! Returns the processor index, in the range [0, ..., nproc-1], of point i.
   !****************************************************************************
   elemental function proc_id_vmu(lo, i)
   
      implicit none
      
      integer :: proc_id_vmu
      type(vmu_layout_type), intent(in) :: lo
      integer, intent(in) :: i

      proc_id_vmu = i / lo%blocksize

   end function proc_id_vmu

   !****************************************************************************
   !                                   IDX_VMU                                  
   !****************************************************************************
   ! Returns the index <idx> in [0, ..., ulim_world] of point (iv,imu,is).
   !****************************************************************************
   elemental function idx_vmu(lo, iv, imu, is)

      implicit none

      integer :: idx_vmu
      type(vmu_layout_type), intent(in) :: lo
      integer, intent(in) :: iv, imu, is

      select case (vms_layout)
      case ('vms')
         idx_vmu = iv - 1 + lo%nvpa * (imu - 1 + lo%nmu * (is - 1))
      case ('mvs')
         idx_vmu = imu - 1 + lo%nmu * (iv - 1 + lo%nvpa * (is - 1))
      end select

   end function idx_vmu

   !****************************************************************************
   !                                 IDX_LOCAL                                  
   !****************************************************************************
   ! When <lo> = <vmu_layout_type> then the interface "idx_local" will call
   ! idx_local_vmu, next "idx_local" will call iz_local_vmu(lo, idx).
   ! Returns true if the point <idx> lies on the current processor <iproc>.
   ! Returns false if the point <idx> does not lie on the current processor <iproc>.
   !****************************************************************************
   elemental function idx_local_vmu(lo, iv, imu, is)

      implicit none
      
      logical :: idx_local_vmu
      type(vmu_layout_type), intent(in) :: lo
      integer, intent(in) :: iv, imu, is

      idx_local_vmu = idx_local(lo, idx(lo, iv, imu, is))
      
   end function idx_local_vmu

   !****************************************************************************
   !                                  IZ_LOCAL                                  
   !****************************************************************************
   ! Returns true if the point <idx=iz> lies on the current processor <iproc>.
   ! Returns false if the point <idx=iz> does not lie on the current processor <iproc>.
   !****************************************************************************
   elemental function iz_local_vmu(lo, iz)
   
      implicit none
      
      logical :: iz_local_vmu
      type(vmu_layout_type), intent(in) :: lo
      integer, intent(in) :: iz

      iz_local_vmu = lo%iproc == proc_id(lo, iz)
      
   end function iz_local_vmu


!###############################################################################
!                                  KYMUS LAYOUT                                 
!###############################################################################
! ....
!###############################################################################

   subroutine init_kymus_layout(nzgrid, ntubes, naky, nakx, nvgrid, nmu, nspec)

      use mp, only: iproc, nproc

      implicit none

      integer, intent(in) :: nzgrid, ntubes, naky, nakx, nvgrid, nmu, nspec
      logical, save :: kymus_initialised = .false.
      
      !-------------------------------------------------------------------------

      ! Only initialise once
      if (kymus_initialised) return
      kymus_initialised = .true.

      ! Parallelisation
      kymus_lo%iproc = iproc
      kymus_lo%nzgrid = nzgrid
      kymus_lo%nzed = 2 * nzgrid + 1
      kymus_lo%ntubes = ntubes
      kymus_lo%naky = naky
      kymus_lo%nakx = nakx
      kymus_lo%nvgrid = nvgrid
      kymus_lo%nvpa = 2 * nvgrid
      kymus_lo%nmu = nmu
      kymus_lo%nspec = nspec
      kymus_lo%llim_world = 0
      kymus_lo%ulim_world = naky * nmu * nspec - 1
      kymus_lo%blocksize = kymus_lo%ulim_world / nproc + 1
      kymus_lo%llim_proc = kymus_lo%blocksize * iproc
      kymus_lo%ulim_proc = min(kymus_lo%ulim_world, kymus_lo%llim_proc + kymus_lo%blocksize - 1)
      kymus_lo%ulim_alloc = max(kymus_lo%llim_proc, kymus_lo%ulim_proc)

   end subroutine init_kymus_layout

   !----------------------------------------------------------------------------
   elemental function is_idx_kymus(lo, i)
      implicit none
      integer :: is_idx_kymus
      type(kymus_layout_type), intent(in) :: lo
      integer, intent(in) :: i
      is_idx_kymus = 1 + mod((i - lo%llim_world) / lo%naky / lo%nmu, lo%nspec)
   end function is_idx_kymus

   !----------------------------------------------------------------------------
   elemental function iky_idx_kymus(lo, i)

      implicit none
      integer :: iky_idx_kymus
      type(kymus_layout_type), intent(in) :: lo
      integer, intent(in) :: i

      select case (kymus_layout)
      case ('kymus')
         iky_idx_kymus = 1 + mod(i - lo%llim_world, lo%naky)
      case ('mukys')
         iky_idx_kymus = 1 + mod((i - lo%llim_world) / lo%nmu, lo%naky)
      end select

   end function iky_idx_kymus

   !----------------------------------------------------------------------------
   elemental function imu_idx_kymus(lo, i)

      implicit none
      integer :: imu_idx_kymus
      type(kymus_layout_type), intent(in) :: lo
      integer, intent(in) :: i

      select case (kymus_layout)
      case ('kymus')
         imu_idx_kymus = 1 + mod((i - lo%llim_world) / lo%naky, lo%nmu)
      case ('mukys')
         imu_idx_kymus = 1 + mod(i - lo%llim_world, lo%nmu)
      end select

   end function imu_idx_kymus

   !----------------------------------------------------------------------------
   elemental function proc_id_kymus(lo, i)
      implicit none
      integer :: proc_id_kymus
      type(kymus_layout_type), intent(in) :: lo
      integer, intent(in) :: i

      proc_id_kymus = i / lo%blocksize

   end function proc_id_kymus

   !----------------------------------------------------------------------------
   elemental function idx_kymus(lo, iky, imu, is)

      implicit none

      integer :: idx_kymus
      type(kymus_layout_type), intent(in) :: lo
      integer, intent(in) :: iky, imu, is

      select case (kymus_layout)
      case ('kymus')
         idx_kymus = iky - 1 + lo%naky * (imu - 1 + lo%nmu * (is - 1))
      case ('mukys')
         idx_kymus = imu - 1 + lo%nmu * (iky - 1 + lo%naky * (is - 1))
      end select

   end function idx_kymus

   !----------------------------------------------------------------------------
   elemental function idx_local_kymus(lo, iky, imu, is)

      implicit none
      logical :: idx_local_kymus
      type(kymus_layout_type), intent(in) :: lo
      integer, intent(in) :: iky, imu, is

      idx_local_kymus = idx_local(lo, idx(lo, iky, imu, is))
   end function idx_local_kymus

   !----------------------------------------------------------------------------
   elemental function iz_local_kymus(lo, iz)
      implicit none
      logical :: iz_local_kymus
      type(kymus_layout_type), intent(in) :: lo
      integer, intent(in) :: iz

      iz_local_kymus = lo%iproc == proc_id(lo, iz)
   end function iz_local_kymus
   

!###############################################################################
!                            SWITCH BETWEEN LAYOUTS                             
!###############################################################################
! The routines "kxkyzidx2vmuidx", "kxyzidx2vmuidx" and "xyzidx2vmuidx" allows us
! to find the correspondence between the points on the different layout types:
!        <ikxkyz> = point (iky,ikx,iz,it,is)
!        <ikxyz>  = point (iy, ikx,iz,it,is)
!        <ixyz>   = point (iy, ix, iz,it,is)
!        <ivmu>   = point (iv,imu,is)
! As well as the corresponding values of (iv,imu,is,iy,ix,iky,ikx,iz,it).
!###############################################################################

   !****************************************************************************
   !                            KXKYZ IDX 2 VMU IDX                             
   !****************************************************************************
   ! From the indices (iv, imu) we can return the point <ivmu> on the vmu_layout.
   ! From the point <ikxkyz> on the kxkyz_layout we can retun (iky, ikx, iz, it).
   !****************************************************************************
   elemental subroutine kxkyzidx2vmuidx(iv, imu, ikxkyz, kxkyz_lo, vmu_lo, iky, ikx, iz, it, ivmu)
   
      implicit none
      
      integer, intent(in) :: iv, imu, ikxkyz
      type(kxkyz_layout_type), intent(in) :: kxkyz_lo
      type(vmu_layout_type), intent(in) :: vmu_lo
      integer, intent(out) :: iky, ikx, iz, it, ivmu

      iky = iky_idx(kxkyz_lo, ikxkyz)
      ikx = ikx_idx(kxkyz_lo, ikxkyz)
      iz = iz_idx(kxkyz_lo, ikxkyz)
      it = it_idx(kxkyz_lo, ikxkyz)
      ivmu = idx(vmu_lo, iv, imu, is_idx(kxkyz_lo, ikxkyz))
      
   end subroutine kxkyzidx2vmuidx

   !****************************************************************************
   !                             KXYZ IDX 2 VMU IDX                             
   !****************************************************************************
   ! From the indices (iv, imu) we can return the point <ivmu> on the vmu_layout.
   ! From the point <ikxyz> on the kxyz_layout we can retun (iy, ikx, iz, it).
   !****************************************************************************
   elemental subroutine kxyzidx2vmuidx(iv, imu, ikxyz, kxyz_lo, vmu_lo, iy, ikx, iz, it, ivmu)
   
      implicit none
      
      integer, intent(in) :: iv, imu, ikxyz
      type(kxyz_layout_type), intent(in) :: kxyz_lo
      type(vmu_layout_type), intent(in) :: vmu_lo
      integer, intent(out) :: iy, ikx, iz, it, ivmu

      iy = iy_idx(kxyz_lo, ikxyz)
      ikx = ikx_idx(kxyz_lo, ikxyz)
      iz = iz_idx(kxyz_lo, ikxyz)
      it = it_idx(kxyz_lo, ikxyz)
      ivmu = idx(vmu_lo, iv, imu, is_idx(kxyz_lo, ikxyz))
      
   end subroutine kxyzidx2vmuidx

   !****************************************************************************
   !                             XYZ IDX 2 VMU IDX                              
   !****************************************************************************
   ! From the indices (iv, imu) we can return the point <ivmu> on the vmu_layout.
   ! From the point <ixyz> on the xyz_layout we can retun (iy, ix, iz, it).
   !****************************************************************************
   elemental subroutine xyzidx2vmuidx(iv, imu, ixyz, xyz_lo, vmu_lo, iy, ix, iz, it, ivmu)
   
      implicit none
      
      integer, intent(in) :: iv, imu, ixyz
      type(xyz_layout_type), intent(in) :: xyz_lo
      type(vmu_layout_type), intent(in) :: vmu_lo
      integer, intent(out) :: iy, ix, iz, it, ivmu

      iy = iy_idx(xyz_lo, ixyz)
      ix = ix_idx(xyz_lo, ixyz)
      iz = iz_idx(xyz_lo, ixyz)
      it = it_idx(xyz_lo, ixyz)
      ivmu = idx(vmu_lo, iv, imu, is_idx(xyz_lo, ixyz))
      
   end subroutine xyzidx2vmuidx

   !****************************************************************************
   !                            KYMUS IDX 2 VMU IDX                             
   !****************************************************************************
   elemental subroutine kymusidx2vmuidx(iv, ikymus, kymus_lo, vmu_lo, iky, ivmu)
   
      implicit none
      
      integer, intent(in) :: iv, ikymus
      type(kymus_layout_type), intent(in) :: kymus_lo
      type(vmu_layout_type), intent(in) :: vmu_lo
      integer, intent(out) :: iky, ivmu

      iky = iky_idx(kymus_lo, ikymus)
      ivmu = idx(vmu_lo, iv, imu_idx(kymus_lo, ikymus), is_idx(kymus_lo, ikymus))
      
   end subroutine kymusidx2vmuidx

!###############################################################################
!                              FINISH LAYOUTS
!###############################################################################

   subroutine finish_layouts

      implicit none

   end subroutine finish_layouts

end module parallelisation_layouts
