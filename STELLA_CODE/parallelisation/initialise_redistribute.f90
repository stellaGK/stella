!###############################################################################
!                                 REDISTRIBUTE                                  
!###############################################################################
! This module initialises the redistribute routines. Before this can be done, the
! parallelisation layouts need to be initialised using init_dist_fn_layouts().
!###############################################################################
module initialise_redistribute

   use redistribute, only: redist_type
   use debug_flags, only: debug => parallelisation_debug

   implicit none

   ! Make routines available to other modules
   public :: init_redistribute, finish_redistribute
   public :: test_kymus_to_vmus_redistribute
   
   ! Make parameters available to other modules
   public :: kxkyz2vmu
   public :: kxyz2vmu
   public :: xyz2vmu
   public :: kymus2vmus

   private

   ! Redistribution information
   type(redist_type) :: kxyz2vmu
   type(redist_type) :: xyz2vmu
   
   ! For the kxkyz_layout, the (mu, vpa) dimensions are local and the (kx, ky, z, tube, s) dimensions are parallelised.
   ! For the vmu_layout, the (kx, ky, z, tube) dimensions are local and the (vpa, mu, s) dimensions are parallelised.
   ! For the mirror term, we need the velocity data to be local (i.e. we need to go from the vmu_layout to the kxkyz_layout)
   ! In other words, we go from  g(naky, nakx, -nzgrid:nzgrid, ntubes, vmu-layout) to gvmu(nvpa, nmu, kxkyz-layout)
   type(redist_type) :: kxkyz2vmu
   
   ! this is a redist_type where we got from a layout with ky, mu and species parallelised
   ! and redistribute so that vpa, mu and species are parallelised
   type(redist_type) :: kymus2vmus
   
   ! Only initialise once
   logical :: initialised_redistribute = .false.
   logical :: initialized_kxkyz_to_vmu_redistribute = .false.
   logical :: initialized_kxyz_to_vmu_redistribute = .false.
   logical :: initialized_xyz_to_vmu_redistribute = .false.
   logical :: initialized_kymus_to_vmus_redistribute = .false.

contains

!###############################################################################
!########################### INITIALISE REDISTRIBUTE ###########################
!###############################################################################

   !****************************************************************************
   !                   Initialise the redistribution routines                   
   !****************************************************************************
   ! Each individual processor needs to know which points of layout A correspond
   ! to the same point on layout B, which is set in this module.
   !****************************************************************************
   subroutine init_redistribute

      use parameters_physics, only: full_flux_surface
      use parameters_physics, only: include_parallel_nonlinearity
      use parameters_numerical, only: split_parallel_dynamics

      implicit none

      !----------------------------------------------------------------------

      ! Only initialise once
      if (initialised_redistribute) return
      initialised_redistribute = .true.

      ! Initialise the redistribute schemes
      call init_kxkyz_to_vmu_redistribute
      if (full_flux_surface) call init_kxyz_to_vmu_redistribute
      if (include_parallel_nonlinearity) call init_xyz_to_vmu_redistribute
      if (.not. split_parallel_dynamics) call init_kymus_to_vmus_redistribute

   end subroutine init_redistribute

   !****************************************************************************
   !                Redistribute from kxkyz-layout to vmu-layout                
   !****************************************************************************
   ! Redistribute the distribution function from gvmu(nvpa, nmu, kxkyz-layout) 
   ! to g(naky, nakx, -nzgrid:nzgrid, ntubes, vmu-layout)
   !****************************************************************************
   subroutine init_kxkyz_to_vmu_redistribute

      use mp, only: nproc
      use parallelisation_layouts, only: kxkyz_lo, vmu_lo
      use parallelisation_layouts, only: iky_idx, ikx_idx, iz_idx, it_idx, is_idx
      use parallelisation_layouts, only: proc_id
      use redistribute, only: index_list_type, init_redist
      use redistribute, only: delete_list, set_redist_character_type
      use grids_velocity, only: nvpa, nmu
      use grids_z, only: nzgrid

      implicit none

      type(index_list_type), dimension(0:nproc - 1) :: to_list, from_list
      integer, dimension(0:nproc - 1) :: nn_to, nn_from
      integer, dimension(3) :: from_low, from_high
      integer, dimension(5) :: to_high, to_low
      integer :: ikxkyz, ip, n, j
      integer :: iv, imu, is, iky, ikx, iz, it
      integer, dimension(:), allocatable :: nlocal
      integer, dimension(:, :), allocatable :: ivmu_local, nn_per_proc
      integer, dimension(:, :, :), allocatable :: proc_of_vmu

      !----------------------------------------------------------------------

      ! Only initialise once
      if (initialized_kxkyz_to_vmu_redistribute) return
      initialized_kxkyz_to_vmu_redistribute = .true.

      ! For each species, find the processor of every (iv, imu, is) point on the
      ! vmu-layout, and the (iv, imu, is) points which are local to this processor
      call get_vmu_points_per_species(kxkyz_lo%nspec, proc_of_vmu, nn_per_proc, nlocal, ivmu_local)

      ! Count number of elements to be redistributed to/from each processor.
      ! The <nn_from> points are the (iv, imu) points of the local ikxkyz-points,
      ! and the <nn_to> points are the local (iv, imu) points of all ikxkyz-points.
      nn_to = 0; nn_from = 0
      do ikxkyz = kxkyz_lo%llim_proc, kxkyz_lo%ulim_proc
         nn_from = nn_from + nn_per_proc(:, is_idx(kxkyz_lo, ikxkyz))
      end do
      do ikxkyz = kxkyz_lo%llim_world, kxkyz_lo%ulim_world
         ip = proc_id(kxkyz_lo, ikxkyz)
         nn_to(ip) = nn_to(ip) + nlocal(is_idx(kxkyz_lo, ikxkyz))
      end do
      
      ! Debug message
      if (debug) then
         write(*,*) ''; write(*,*) 'Initialise kxkyz_to_vmu_redistribute:'
         do ip = 0, nproc - 1
            write(*,'(A,I0,A,I0,A,I0)') '    kxkyz_to_vmu on proc ', ip, ': nn_from = ', nn_from(ip), ', nn_to = ', nn_to(ip)
         end do
      end if

      ! We go from (nvpa, nmu, kxkyz-layout) to (nky, nkx, nz, ntubes, vmu-layout)
      ! For each of the <nn_from> points which are local to the processor, we need
      ! to tell it what (ivpa, imu, ikxkyz) point it is, similarly, for each <nn_to> 
      ! point, we need to tell it what (iky, ikx, iz, itube, ivmu) point it is.
      ! Here <from_list> and <to_list> will hold the indices of these points.
      do ip = 0, nproc - 1
         if (nn_from(ip) > 0) then
            allocate (from_list(ip)%first(nn_from(ip)))     ! nvpa
            allocate (from_list(ip)%second(nn_from(ip)))    ! nmu
            allocate (from_list(ip)%third(nn_from(ip)))     ! kxkyz-layout
         end if
         if (nn_to(ip) > 0) then
            allocate (to_list(ip)%first(nn_to(ip)))         ! nky
            allocate (to_list(ip)%second(nn_to(ip)))        ! nkx
            allocate (to_list(ip)%third(nn_to(ip)))         ! nz
            allocate (to_list(ip)%fourth(nn_to(ip)))        ! ntubes
            allocate (to_list(ip)%fifth(nn_to(ip)))         ! vmu-layout
         end if
      end do

      ! Get local indices of elements distributed to/from other processors.
      ! Both lists are ordered as a loop over (ikxkyz, imu, iv) with iv the fastest
      ! index, so that the order of the points agrees between the sending and receiving processor.
      nn_to = 0
      nn_from = 0

      ! The ikxkyz-points local to this processor are the <nn_from> points,
      ! for which we set the (vpa, mu, ikxkyz) indices
      do ikxkyz = kxkyz_lo%llim_proc, kxkyz_lo%ulim_proc
         is = is_idx(kxkyz_lo, ikxkyz)
         do imu = 1, nmu
            do iv = 1, nvpa
               ip = proc_of_vmu(iv, imu, is)
               n = nn_from(ip) + 1
               nn_from(ip) = n
               from_list(ip)%first(n) = iv
               from_list(ip)%second(n) = imu
               from_list(ip)%third(n) = ikxkyz
            end do
         end do
      end do

      ! The ivmu-points local to this processor are the <nn_to> points,
      ! for which we set the (ky, kx, z, tube, ivmu) indices
      do ikxkyz = kxkyz_lo%llim_world, kxkyz_lo%ulim_world
         is = is_idx(kxkyz_lo, ikxkyz)
         if (nlocal(is) == 0) cycle
         ip = proc_id(kxkyz_lo, ikxkyz)
         iky = iky_idx(kxkyz_lo, ikxkyz)
         ikx = ikx_idx(kxkyz_lo, ikxkyz)
         iz = iz_idx(kxkyz_lo, ikxkyz)
         it = it_idx(kxkyz_lo, ikxkyz)
         n = nn_to(ip)
         do j = 1, nlocal(is)
            to_list(ip)%first(n + j) = iky
            to_list(ip)%second(n + j) = ikx
            to_list(ip)%third(n + j) = iz
            to_list(ip)%fourth(n + j) = it
            to_list(ip)%fifth(n + j) = ivmu_local(j, is)
         end do
         nn_to(ip) = n + nlocal(is)
      end do
      
      deallocate (proc_of_vmu, nn_per_proc, nlocal, ivmu_local)

      ! Redistribute the distribution function from gvmu(nvpa, nmu, kxkyz-layout) 
      ! Here, we set the lower and upper limits of the dimensions
      from_low(1) = 1
      from_low(2) = 1
      from_low(3) = kxkyz_lo%llim_proc
      
      from_high(1) = nvpa
      from_high(2) = nmu
      from_high(3) = kxkyz_lo%ulim_alloc

      ! Redistribute the distribution function to g(naky, nakx, -nzgrid:nzgrid, ntubes, vmu-layout)
      ! Here, we set the lower and upper limits of the dimensions
      to_low(1) = 1
      to_low(2) = 1
      to_low(3) = -nzgrid
      to_low(4) = 1
      to_low(5) = vmu_lo%llim_proc
      
      to_high(1) = vmu_lo%naky
      to_high(2) = vmu_lo%nakx
      to_high(3) = vmu_lo%nzgrid
      to_high(4) = vmu_lo%ntubes
      to_high(5) = vmu_lo%ulim_alloc

      ! Send <to_list>, <from_list> and its limits to each processor
      call set_redist_character_type(kxkyz2vmu, 'kxkyz2vmu')
      call init_redist(kxkyz2vmu, 'c', to_low, to_high, to_list, from_low, from_high, from_list)
      call delete_list(to_list)
      call delete_list(from_list)

   end subroutine init_kxkyz_to_vmu_redistribute

   !****************************************************************************
   !                                      Title
   !****************************************************************************
   subroutine init_kxyz_to_vmu_redistribute

      use mp, only: nproc
      use parallelisation_layouts, only: kxyz_lo, vmu_lo
      use parallelisation_layouts, only: iy_idx, ikx_idx, iz_idx, it_idx, is_idx
      use parallelisation_layouts, only: proc_id
      use redistribute, only: index_list_type, init_redist
      use redistribute, only: delete_list, set_redist_character_type
      use grids_velocity, only: nvpa, nmu
      use grids_z, only: nzgrid

      implicit none

      type(index_list_type), dimension(0:nproc - 1) :: to_list, from_list
      integer, dimension(0:nproc - 1) :: nn_to, nn_from
      integer, dimension(3) :: from_low, from_high
      integer, dimension(5) :: to_high, to_low
      integer :: ikxyz
      integer :: iv, imu, is, iy, ikx, iz, it
      integer :: ip, n, j
      integer, dimension(:), allocatable :: nlocal
      integer, dimension(:, :), allocatable :: ivmu_local, nn_per_proc
      integer, dimension(:, :, :), allocatable :: proc_of_vmu

      !----------------------------------------------------------------------

      ! Only initialise once
      if (initialized_kxyz_to_vmu_redistribute) return
      initialized_kxyz_to_vmu_redistribute = .true.

      ! For each species, find the processor of every (iv, imu, is) point on the
      ! vmu-layout, and the (iv, imu, is) points which are local to this processor
      call get_vmu_points_per_species(kxyz_lo%nspec, proc_of_vmu, nn_per_proc, nlocal, ivmu_local)

      ! Count number of elements to be redistributed to/from each processor
      nn_to = 0
      nn_from = 0
      do ikxyz = kxyz_lo%llim_proc, kxyz_lo%ulim_proc
         nn_from = nn_from + nn_per_proc(:, is_idx(kxyz_lo, ikxyz))
      end do
      do ikxyz = kxyz_lo%llim_world, kxyz_lo%ulim_world
         ip = proc_id(kxyz_lo, ikxyz)
         nn_to(ip) = nn_to(ip) + nlocal(is_idx(kxyz_lo, ikxyz))
      end do

      do ip = 0, nproc - 1
         if (nn_from(ip) > 0) then
            allocate (from_list(ip)%first(nn_from(ip)))
            allocate (from_list(ip)%second(nn_from(ip)))
            allocate (from_list(ip)%third(nn_from(ip)))
         end if
         if (nn_to(ip) > 0) then
            allocate (to_list(ip)%first(nn_to(ip)))
            allocate (to_list(ip)%second(nn_to(ip)))
            allocate (to_list(ip)%third(nn_to(ip)))
            allocate (to_list(ip)%fourth(nn_to(ip)))
            allocate (to_list(ip)%fifth(nn_to(ip)))
         end if
      end do

      ! Get local indices of elements distributed to/from other processors.
      ! Both lists are ordered as a loop over (ikxyz, imu, iv) with iv the fastest
      ! index, so that the order of the points agrees between the sending and receiving processor.
      nn_to = 0
      nn_from = 0

      ! if ikxyz index local, set:
      ! ip = corresponding vmu processor
      ! from_list%first-third arrays = iv,imu,ikxyz  (ie vmu indices)
      ! later will send from_list to proc ip
      do ikxyz = kxyz_lo%llim_proc, kxyz_lo%ulim_proc
         is = is_idx(kxyz_lo, ikxyz)
         do imu = 1, nmu
            do iv = 1, nvpa
               ip = proc_of_vmu(iv, imu, is)
               n = nn_from(ip) + 1
               nn_from(ip) = n
               from_list(ip)%first(n) = iv
               from_list(ip)%second(n) = imu
               from_list(ip)%third(n) = ikxyz
            end do
         end do
      end do

      ! if vmu index local, set ip to corresponding y processor
      ! set to_list%first-fifth arrays = iy,ikx,iz,it,ivmu  (ie y indices)
      ! will receive to_list from ip
      do ikxyz = kxyz_lo%llim_world, kxyz_lo%ulim_world
         is = is_idx(kxyz_lo, ikxyz)
         if (nlocal(is) == 0) cycle
         ip = proc_id(kxyz_lo, ikxyz)
         iy = iy_idx(kxyz_lo, ikxyz)
         ikx = ikx_idx(kxyz_lo, ikxyz)
         iz = iz_idx(kxyz_lo, ikxyz)
         it = it_idx(kxyz_lo, ikxyz)
         n = nn_to(ip)
         do j = 1, nlocal(is)
            to_list(ip)%first(n + j) = iy
            to_list(ip)%second(n + j) = ikx
            to_list(ip)%third(n + j) = iz
            to_list(ip)%fourth(n + j) = it
            to_list(ip)%fifth(n + j) = ivmu_local(j, is)
         end do
         nn_to(ip) = n + nlocal(is)
      end do

      deallocate (proc_of_vmu, nn_per_proc, nlocal, ivmu_local)

      from_low(1) = 1
      from_low(2) = 1
      from_low(3) = kxyz_lo%llim_proc

      from_high(1) = nvpa
      from_high(2) = nmu
      from_high(3) = kxyz_lo%ulim_alloc

      to_low(1) = 1
      to_low(2) = 1
      to_low(3) = -nzgrid
      to_low(4) = 1
      to_low(5) = vmu_lo%llim_proc

      to_high(1) = vmu_lo%ny
      to_high(2) = vmu_lo%nakx / 2 + 1
      to_high(3) = vmu_lo%nzed
      to_high(4) = vmu_lo%ntubes
      to_high(5) = vmu_lo%ulim_alloc

      call set_redist_character_type(kxyz2vmu, 'kxyz2vmu')
      call init_redist(kxyz2vmu, 'c', to_low, to_high, to_list, from_low, from_high, from_list)
      call delete_list(to_list)
      call delete_list(from_list)

   end subroutine init_kxyz_to_vmu_redistribute

   !****************************************************************************
   !                                      Title
   !****************************************************************************
   subroutine init_xyz_to_vmu_redistribute

      use mp, only: nproc
      use parallelisation_layouts, only: xyz_lo, vmu_lo
      use parallelisation_layouts, only: iy_idx, ix_idx, iz_idx, it_idx, is_idx
      use parallelisation_layouts, only: proc_id
      use redistribute, only: index_list_type, init_redist
      use redistribute, only: delete_list, set_redist_character_type
      use grids_velocity, only: nvpa, nmu
      use grids_z, only: nzgrid

      implicit none

      type(index_list_type), dimension(0:nproc - 1) :: to_list, from_list
      integer, dimension(0:nproc - 1) :: nn_to, nn_from
      integer, dimension(3) :: from_low, from_high
      integer, dimension(5) :: to_high, to_low
      integer :: ixyz
      integer :: iv, imu, is, iy, ix, iz, it
      integer :: ip, n, j
      integer, dimension(:), allocatable :: nlocal
      integer, dimension(:, :), allocatable :: ivmu_local, nn_per_proc
      integer, dimension(:, :, :), allocatable :: proc_of_vmu

      !----------------------------------------------------------------------

      ! Only initialise once
      if (initialized_xyz_to_vmu_redistribute) return
      initialized_xyz_to_vmu_redistribute = .true.

      ! For each species, find the processor of every (iv, imu, is) point on the
      ! vmu-layout, and the (iv, imu, is) points which are local to this processor
      call get_vmu_points_per_species(xyz_lo%nspec, proc_of_vmu, nn_per_proc, nlocal, ivmu_local)

      ! Count number of elements to be redistributed to/from each processor
      nn_to = 0
      nn_from = 0
      do ixyz = xyz_lo%llim_proc, xyz_lo%ulim_proc
         nn_from = nn_from + nn_per_proc(:, is_idx(xyz_lo, ixyz))
      end do
      do ixyz = xyz_lo%llim_world, xyz_lo%ulim_world
         ip = proc_id(xyz_lo, ixyz)
         nn_to(ip) = nn_to(ip) + nlocal(is_idx(xyz_lo, ixyz))
      end do

      do ip = 0, nproc - 1
         if (nn_from(ip) > 0) then
            allocate (from_list(ip)%first(nn_from(ip)))
            allocate (from_list(ip)%second(nn_from(ip)))
            allocate (from_list(ip)%third(nn_from(ip)))
         end if
         if (nn_to(ip) > 0) then
            allocate (to_list(ip)%first(nn_to(ip)))
            allocate (to_list(ip)%second(nn_to(ip)))
            allocate (to_list(ip)%third(nn_to(ip)))
            allocate (to_list(ip)%fourth(nn_to(ip)))
            allocate (to_list(ip)%fifth(nn_to(ip)))
         end if
      end do

      ! Get local indices of elements distributed to/from other processors.
      ! Both lists are ordered as a loop over (ixyz, imu, iv) with iv the fastest
      ! index, so that the order of the points agrees between the sending and receiving processor.
      nn_to = 0
      nn_from = 0

      ! if ixyz index local, set:
      ! ip = corresponding vmu processor
      ! from_list%first-third arrays = iv,imu,ixyz  (ie vmu indices)
      ! later will send from_list to proc ip
      do ixyz = xyz_lo%llim_proc, xyz_lo%ulim_proc
         is = is_idx(xyz_lo, ixyz)
         do imu = 1, nmu
            do iv = 1, nvpa
               ip = proc_of_vmu(iv, imu, is)
               n = nn_from(ip) + 1
               nn_from(ip) = n
               from_list(ip)%first(n) = iv
               from_list(ip)%second(n) = imu
               from_list(ip)%third(n) = ixyz
            end do
         end do
      end do

      ! if vmu index local, set ip to corresponding y processor
      ! set to_list%first-fifth arrays = iy,ix,iz,it,ivmu  (ie y indices)
      ! will receive to_list from ip
      do ixyz = xyz_lo%llim_world, xyz_lo%ulim_world
         is = is_idx(xyz_lo, ixyz)
         if (nlocal(is) == 0) cycle
         ip = proc_id(xyz_lo, ixyz)
         iy = iy_idx(xyz_lo, ixyz)
         ix = ix_idx(xyz_lo, ixyz)
         iz = iz_idx(xyz_lo, ixyz)
         it = it_idx(xyz_lo, ixyz)
         n = nn_to(ip)
         do j = 1, nlocal(is)
            to_list(ip)%first(n + j) = iy
            to_list(ip)%second(n + j) = ix
            to_list(ip)%third(n + j) = iz
            to_list(ip)%fourth(n + j) = it
            to_list(ip)%fifth(n + j) = ivmu_local(j, is)
         end do
         nn_to(ip) = n + nlocal(is)
      end do

      deallocate (proc_of_vmu, nn_per_proc, nlocal, ivmu_local)

      from_low(1) = 1
      from_low(2) = 1
      from_low(3) = xyz_lo%llim_proc

      from_high(1) = nvpa
      from_high(2) = nmu
      from_high(3) = xyz_lo%ulim_alloc

      to_low(1) = 1
      to_low(2) = 1
      to_low(3) = -nzgrid
      to_low(4) = 1
      to_low(5) = vmu_lo%llim_proc

      to_high(1) = vmu_lo%ny
      to_high(2) = vmu_lo%nx
      to_high(3) = vmu_lo%nzed
      to_high(4) = vmu_lo%ntubes
      to_high(5) = vmu_lo%ulim_alloc

      call set_redist_character_type(xyz2vmu, 'xyz2vmu')
      call init_redist(xyz2vmu, 'r', to_low, to_high, to_list, from_low, from_high, from_list)
      call delete_list(to_list)
      call delete_list(from_list)

   end subroutine init_xyz_to_vmu_redistribute

   !****************************************************************************
   !                                      Title
   !****************************************************************************
   subroutine init_kymus_to_vmus_redistribute

      use mp, only: nproc
      use parallelisation_layouts, only: kymus_lo, vmu_lo
      use parallelisation_layouts, only: kymusidx2vmuidx
      use parallelisation_layouts, only: idx_local, proc_id
      use redistribute, only: index_list_type, init_redist
      use redistribute, only: delete_list, set_redist_character_type
      use grids_velocity, only: nvpa
      use grids_kxky, only: nakx
      use grids_z, only: nzgrid, ntubes

      implicit none

      type(index_list_type), dimension(0:nproc - 1) :: to_list, from_list
      integer, dimension(0:nproc - 1) :: nn_to, nn_from
      integer, dimension(5) :: from_low, from_high
      integer, dimension(5) :: to_high, to_low
      integer :: ikymus, ivmu
      integer :: iv, iky, ikx, iz, it
      integer :: ip, n, nkxzt

      !----------------------------------------------------------------------

      ! Only initialise once
      if (initialized_kymus_to_vmus_redistribute) return
      initialized_kymus_to_vmus_redistribute = .true.

      ! The (ikx, iz, it) points do not change the processor of the point, so each
      ! (ikymus, iv) point corresponds to <nkxzt> points that need to be redistributed
      nkxzt = nakx * (2 * nzgrid + 1) * ntubes

      ! Count number of elements to be redistributed to/from each processor
      nn_to = 0
      nn_from = 0
      do ikymus = kymus_lo%llim_world, kymus_lo%ulim_world
         do iv = 1, nvpa
            call kymusidx2vmuidx(iv, ikymus, kymus_lo, vmu_lo, iky, ivmu)
            if (idx_local(kymus_lo, ikymus)) &
                 nn_from(proc_id(vmu_lo, ivmu)) = nn_from(proc_id(vmu_lo, ivmu)) + nkxzt
            if (idx_local(vmu_lo, ivmu)) &
                 nn_to(proc_id(kymus_lo, ikymus)) = nn_to(proc_id(kymus_lo, ikymus)) + nkxzt
         end do
      end do
            
      do ip = 0, nproc - 1
         if (nn_from(ip) > 0) then
            allocate (from_list(ip)%first(nn_from(ip)))
            allocate (from_list(ip)%second(nn_from(ip)))
            allocate (from_list(ip)%third(nn_from(ip)))
            allocate (from_list(ip)%fourth(nn_from(ip)))
            allocate (from_list(ip)%fifth(nn_from(ip)))
         end if
         if (nn_to(ip) > 0) then
            allocate (to_list(ip)%first(nn_to(ip)))
            allocate (to_list(ip)%second(nn_to(ip)))
            allocate (to_list(ip)%third(nn_to(ip)))
            allocate (to_list(ip)%fourth(nn_to(ip)))
            allocate (to_list(ip)%fifth(nn_to(ip)))
         end if
      end do

      ! Get local indices of elements distributed to/from other processors
      nn_to = 0
      nn_from = 0

      ! Loop over all indices in the kymus layout and find the corresponding indices for iky and ivmu
      ! Only the (ikymus, iv) points with ikymus or ivmu local to this processor need the inner loops
      do ikymus = kymus_lo%llim_world, kymus_lo%ulim_world
         do iv = 1, nvpa
         
            ! Obtain corresponding ky indices
            call kymusidx2vmuidx(iv, ikymus, kymus_lo, vmu_lo, iky, ivmu)
            
            ! if kymus index local, set:
            ! ip = corresponding y processor
            ! from_list%first-third arrays = ikx,iz,it,iv,ikymus
            ! later will send from_list to proc ip
            if (idx_local(kymus_lo, ikymus)) then
               ip = proc_id(vmu_lo, ivmu)
               n = nn_from(ip)
               do it = 1, ntubes
                  do iz = -nzgrid, nzgrid
                     do ikx = 1, nakx
                        n = n + 1
                        from_list(ip)%first(n) = ikx
                        from_list(ip)%second(n) = iz
                        from_list(ip)%third(n) = it
                        from_list(ip)%fourth(n) = iv
                        from_list(ip)%fifth(n) = ikymus
                     end do
                  end do
               end do
               nn_from(ip) = n
            end if
            if (idx_local(vmu_lo, ivmu)) then
               ip = proc_id(kymus_lo, ikymus)
               n = nn_to(ip)
               do it = 1, ntubes
                  do iz = -nzgrid, nzgrid
                     do ikx = 1, nakx
                        n = n + 1
                        to_list(ip)%first(n) = iky
                        to_list(ip)%second(n) = ikx
                        to_list(ip)%third(n) = iz
                        to_list(ip)%fourth(n) = it
                        to_list(ip)%fifth(n) = ivmu
                     end do
                  end do
               end do
               nn_to(ip) = n
            end if
            
         end do
      end do

      from_low(1) = 1
      from_low(2) = -kymus_lo%nzgrid
      from_low(3) = 1
      from_low(4) = 1
      from_low(5) = kymus_lo%llim_proc

      from_high(1) = kymus_lo%nakx
      from_high(2) = kymus_lo%nzgrid
      from_high(3) = kymus_lo%ntubes
      from_high(4) = kymus_lo%nvpa
      from_high(5) = kymus_lo%ulim_alloc

      to_low(1) = 1
      to_low(2) = 1
      to_low(3) = -vmu_lo%nzgrid
      to_low(4) = 1
      to_low(5) = vmu_lo%llim_proc

      to_high(1) = vmu_lo%naky
      to_high(2) = vmu_lo%nakx
      to_high(3) = vmu_lo%nzgrid
      to_high(4) = vmu_lo%ntubes
      to_high(5) = vmu_lo%ulim_alloc

      call set_redist_character_type(kymus2vmus, 'kymus2vmus')
      call init_redist(kymus2vmus, 'c', to_low, to_high, to_list, from_low, from_high, from_list)
      call delete_list(to_list)
      call delete_list(from_list)

   end subroutine init_kymus_to_vmus_redistribute

   
!###############################################################################
!################################# UTILITIES ###################################
!###############################################################################

   !****************************************************************************
   !                    Points on the vmu-layout per species                    
   !****************************************************************************
   ! The (iv, imu) dimensions are local on the kxkyz, kxyz and xyz layouts, while
   ! the species index is shared with the vmu-layout. Therefore, the mapping from
   ! (iv, imu, is) to the vmu-layout only needs to be computed once per species,
   ! instead of for every point of the full (kx, ky, z, tube, s, mu, vpa) grid.
   ! 
   ! <proc_of_vmu(iv, imu, is)> is the processor of the point (iv, imu, is)
   ! <nn_per_proc(ip, is)> is the number of (iv, imu) points of species <is> on processor <ip>
   ! <ivmu_local(1:nlocal(is), is)> are the ivmu-points of species <is> on this processor,
   ! ordered as a loop over (imu, iv) with iv the fastest index.
   !****************************************************************************
   subroutine get_vmu_points_per_species(nspec, proc_of_vmu, nn_per_proc, nlocal, ivmu_local)

      use mp, only: nproc
      use parallelisation_layouts, only: vmu_lo
      use parallelisation_layouts, only: idx, idx_local, proc_id
      use grids_velocity, only: nvpa, nmu

      implicit none

      integer, intent(in) :: nspec
      integer, dimension(:, :, :), allocatable, intent(out) :: proc_of_vmu
      integer, dimension(:, :), allocatable, intent(out) :: nn_per_proc, ivmu_local
      integer, dimension(:), allocatable, intent(out) :: nlocal
      integer :: iv, imu, is, ivmu, ip

      !----------------------------------------------------------------------

      allocate (proc_of_vmu(nvpa, nmu, nspec))
      allocate (nn_per_proc(0:nproc - 1, nspec))
      allocate (ivmu_local(nvpa * nmu, nspec))
      allocate (nlocal(nspec))

      nn_per_proc = 0
      nlocal = 0
      do is = 1, nspec
         do imu = 1, nmu
            do iv = 1, nvpa
               ivmu = idx(vmu_lo, iv, imu, is)
               ip = proc_id(vmu_lo, ivmu)
               proc_of_vmu(iv, imu, is) = ip
               nn_per_proc(ip, is) = nn_per_proc(ip, is) + 1
               if (idx_local(vmu_lo, ivmu)) then
                  nlocal(is) = nlocal(is) + 1
                  ivmu_local(nlocal(is), is) = ivmu
               end if
            end do
         end do
      end do

   end subroutine get_vmu_points_per_species

!###############################################################################
!############################ FINISH REDISTRIBUTE ##############################
!###############################################################################

   subroutine finish_redistribute

      implicit none

      initialised_redistribute = .false.

   end subroutine finish_redistribute
   
!###############################################################################
!############################## TESTS - DEBUGGING ##############################
!###############################################################################

   subroutine test_kymus_to_vmus_redistribute

      use redistribute, only: scatter, gather
      use arrays_distribution_function, only: g_kymus, gnew
      use grids_velocity, only: nvpa
      use mp, only: proc0, send, receive
      use grids_z, only: nzgrid, ntubes
      use grids_kxky, only: nakx, naky
      use parallelisation_layouts, only: vmu_lo, kymus_lo
      use parallelisation_layouts, only: iv_idx, imu_idx, is_idx, iky_idx, idx_local, proc_id

      implicit none

      complex, dimension (:, :, :, :), allocatable :: gtmp
      integer :: ikymus, ivmu, iv, imu, is, it, iz, ikx, iky, izmod

      !----------------------------------------------------------------------

      ! Allocate temporary arrays
      allocate (gtmp(naky, nakx, 2*nzgrid+1, ntubes))
      
      do ivmu = vmu_lo%llim_world, vmu_lo%ulim_world
         iv = iv_idx(vmu_lo, ivmu)
         imu = imu_idx(vmu_lo, ivmu)
         is = is_idx(vmu_lo, ivmu)
         
         if (idx_local(vmu_lo, iv, imu, is)) then
            do iz = -nzgrid, nzgrid
               izmod = iz + nzgrid + 1
               gtmp(:, :, izmod, :) = gnew(:, :, iz, :, ivmu)
            end do
            if (.not. proc0) then
               call send (gtmp, 0)
            end if
         else if (proc0) then
            call receive (gtmp, proc_id(vmu_lo, ivmu))
         end if
         
         if (proc0) then
            do it = 1, ntubes
               do iz = 1, 2*nzgrid+1
                  do ikx = 1, nakx
                     do iky = 1, naky
                        write (26, *) ikx, iky, iz, it, iv, imu, is, real(gtmp(iky, ikx, iz, it)), aimag(gtmp(iky, ikx, iz, it))
                     end do
                  end do
               end do
            end do
         end if
         
      end do
      
      ! Deallocate temporary arrays
      deallocate (gtmp)

      ! Test scatter
      call scatter(kymus2vmus, gnew, g_kymus)
      
      ! Allocate temporary arrays
      allocate (gtmp(nakx, 2*nzgrid+1, ntubes, nvpa))
      do ikymus = kymus_lo%llim_world, kymus_lo%ulim_world
         iky = iky_idx(kymus_lo, ikymus)
         imu = imu_idx(kymus_lo, ikymus)
         is = is_idx(kymus_lo, ikymus)
         if (idx_local(kymus_lo, iky, imu, is)) then
            do iz = -nzgrid, nzgrid
               izmod = iz + nzgrid + 1
               gtmp(:, izmod, :, :) = g_kymus(:, iz, :, :, ivmu)
            end do
            if (.not. proc0) then
               call send (gtmp, 0)
            end if
         else if (proc0) then
            call receive (gtmp, proc_id(kymus_lo, ikymus))
         end if
         if (proc0) then
            do iv = 1, nvpa
               do it = 1, ntubes
                  do iz = 1, 2*nzgrid+1
                     do ikx = 1, nakx
                        write (27, *) ikx, iky, iz, it, iv, imu, is, real(gtmp(ikx, iz, it, iv)), aimag(gtmp(ikx, iz, it, iv))
                     end do
                  end do
               end do
            end do
         end if
      end do
      
      ! Deallocate temporary arrays
      deallocate (gtmp)

      ! Test Gather
      call gather(kymus2vmus, g_kymus, gnew)
      
      ! Allocate temporary arrays
      allocate (gtmp(naky, nakx, 2*nzgrid+1, ntubes))
      
      do ivmu = vmu_lo%llim_world, vmu_lo%ulim_world
         iv = iv_idx(vmu_lo, ivmu)
         imu = imu_idx(vmu_lo, ivmu)
         is = is_idx(vmu_lo, ivmu)
         if (idx_local(vmu_lo, iv, imu, is)) then
            do iz = -nzgrid, nzgrid
               izmod = iz + nzgrid + 1
               gtmp(:, :, izmod, :) = gnew(:, :, iz, :, ivmu)
            end do
            if (.not. proc0) then
               call send (gtmp, 0)
            end if
         else if (proc0) then
            call receive (gtmp, proc_id(vmu_lo, ivmu))
         end if
         if (proc0) then
            do it = 1, ntubes
               do iz = 1, 2*nzgrid+1
                  do ikx = 1, nakx
                     do iky = 1, naky
                        write (28, *) ikx, iky, iz, it, iv, imu, is, real(gtmp(iky, ikx, iz, it)), aimag(gtmp(iky, ikx, iz, it))
                     end do
                  end do
               end do
            end do
         end if
      end do
      
      ! Deallocate temporary arrays
      deallocate (gtmp)
     
   end subroutine test_kymus_to_vmus_redistribute

end module initialise_redistribute
