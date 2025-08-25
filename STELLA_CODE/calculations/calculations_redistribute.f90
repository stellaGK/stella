module calculations_redistribute

   use redistribute, only: redist_type

   implicit none

   public :: init_redistribute, finish_redistribute
   public :: test_kymus_to_vmus_redistribute
   public :: kxkyz2vmu
   public :: kxyz2vmu
   public :: xyz2vmu
   public :: kymus2vmus

   private

   type(redist_type) :: kxkyz2vmu
   type(redist_type) :: kxyz2vmu
   type(redist_type) :: xyz2vmu
   ! this is a redist_type where we got from a layout with ky, mu and species parallelised
   ! and redistribute so that vpa, mu and species are parallelised
   type(redist_type) :: kymus2vmus
   
   logical :: redistribute_initialized = .false.

contains

   subroutine init_redistribute

      use parameters_physics, only: full_flux_surface
      use parameters_physics, only: include_parallel_nonlinearity
      use parameters_numerical, only: split_parallel_dynamics

      implicit none

      if (redistribute_initialized) return
      redistribute_initialized = .true.

      call init_kxkyz_to_vmu_redistribute
      if (full_flux_surface) call init_kxyz_to_vmu_redistribute
      if (include_parallel_nonlinearity) call init_xyz_to_vmu_redistribute
      if (.not.split_parallel_dynamics) call init_kymus_to_vmus_redistribute

   end subroutine init_redistribute

   subroutine init_kxkyz_to_vmu_redistribute

      use mp, only: nproc
      use stella_layouts, only: kxkyz_lo, vmu_lo
      use stella_layouts, only: kxkyzidx2vmuidx
      use stella_layouts, only: idx_local, proc_id
      use redistribute, only: index_list_type, init_redist
      use redistribute, only: delete_list, set_redist_character_type
      use grids_velocity, only: nvpa, nmu
      use grids_z, only: nzgrid

      implicit none

      type(index_list_type), dimension(0:nproc - 1) :: to_list, from_list
      integer, dimension(0:nproc - 1) :: nn_to, nn_from
      integer, dimension(3) :: from_low, from_high
      integer, dimension(5) :: to_high, to_low
      integer :: ikxkyz, ivmu
      integer :: iv, imu, iky, ikx, iz, it
      integer :: ip, n
      logical :: initialized = .false.

      if (initialized) return
      initialized = .true.

      ! count number of elements to be redistributed to/from each processor
      nn_to = 0
      nn_from = 0
      do ikxkyz = kxkyz_lo%llim_world, kxkyz_lo%ulim_world
         do imu = 1, nmu
            do iv = 1, nvpa
               call kxkyzidx2vmuidx(iv, imu, ikxkyz, kxkyz_lo, vmu_lo, iky, ikx, iz, it, ivmu)
               if (idx_local(kxkyz_lo, ikxkyz)) &
                  nn_from(proc_id(vmu_lo, ivmu)) = nn_from(proc_id(vmu_lo, ivmu)) + 1
               if (idx_local(vmu_lo, ivmu)) &
                  nn_to(proc_id(kxkyz_lo, ikxkyz)) = nn_to(proc_id(kxkyz_lo, ikxkyz)) + 1
            end do
         end do
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

      ! get local indices of elements distributed to/from other processors
      nn_to = 0
      nn_from = 0

      ! loop over all vmu indices, find corresponding y indices
      do ikxkyz = kxkyz_lo%llim_world, kxkyz_lo%ulim_world
         do imu = 1, nmu
            do iv = 1, nvpa
               ! obtain corresponding y indices
               call kxkyzidx2vmuidx(iv, imu, ikxkyz, kxkyz_lo, vmu_lo, iky, ikx, iz, it, ivmu)
               ! if vmu index local, set:
               ! ip = corresponding y processor
               ! from_list%first-third arrays = iv,imu,ikxkyz  (ie vmu indices)
               ! later will send from_list to proc ip
               if (idx_local(kxkyz_lo, ikxkyz)) then
                  ip = proc_id(vmu_lo, ivmu)
                  n = nn_from(ip) + 1
                  nn_from(ip) = n
                  from_list(ip)%first(n) = iv
                  from_list(ip)%second(n) = imu
                  from_list(ip)%third(n) = ikxkyz
               end if
               ! if y index local, set ip to corresponding vmu processor
               ! set to_list%first,second arrays = iky,iy  (ie y indices)
               ! will receive to_list from ip
               if (idx_local(vmu_lo, ivmu)) then
                  ip = proc_id(kxkyz_lo, ikxkyz)
                  n = nn_to(ip) + 1
                  nn_to(ip) = n
                  to_list(ip)%first(n) = iky
                  to_list(ip)%second(n) = ikx
                  to_list(ip)%third(n) = iz
                  to_list(ip)%fourth(n) = it
                  to_list(ip)%fifth(n) = ivmu
               end if
            end do
         end do
      end do

      from_low(1) = 1
      from_low(2) = 1
      from_low(3) = kxkyz_lo%llim_proc

      from_high(1) = nvpa
      from_high(2) = nmu
      from_high(3) = kxkyz_lo%ulim_alloc

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

      call set_redist_character_type(kxkyz2vmu, 'kxkyz2vmu')

      call init_redist(kxkyz2vmu, 'c', to_low, to_high, to_list, &
                       from_low, from_high, from_list)

      call delete_list(to_list)
      call delete_list(from_list)

   end subroutine init_kxkyz_to_vmu_redistribute

   subroutine init_kxyz_to_vmu_redistribute

      use mp, only: nproc
      use stella_layouts, only: kxyz_lo, vmu_lo
      use stella_layouts, only: kxyzidx2vmuidx
      use stella_layouts, only: idx_local, proc_id
      use redistribute, only: index_list_type, init_redist
      use redistribute, only: delete_list, set_redist_character_type
      use grids_velocity, only: nvpa, nmu
      use grids_z, only: nzgrid

      implicit none

      type(index_list_type), dimension(0:nproc - 1) :: to_list, from_list
      integer, dimension(0:nproc - 1) :: nn_to, nn_from
      integer, dimension(3) :: from_low, from_high
      integer, dimension(5) :: to_high, to_low
      integer :: ikxyz, ivmu
      integer :: iv, imu, iy, ikx, iz, it
      integer :: ip, n
      logical :: initialized = .false.

      if (initialized) return
      initialized = .true.

      ! count number of elements to be redistributed to/from each processor
      nn_to = 0
      nn_from = 0
      do ikxyz = kxyz_lo%llim_world, kxyz_lo%ulim_world
         do imu = 1, nmu
            do iv = 1, nvpa
               call kxyzidx2vmuidx(iv, imu, ikxyz, kxyz_lo, vmu_lo, iy, ikx, iz, it, ivmu)
               if (idx_local(kxyz_lo, ikxyz)) &
                  nn_from(proc_id(vmu_lo, ivmu)) = nn_from(proc_id(vmu_lo, ivmu)) + 1
               if (idx_local(vmu_lo, ivmu)) &
                  nn_to(proc_id(kxyz_lo, ikxyz)) = nn_to(proc_id(kxyz_lo, ikxyz)) + 1
            end do
         end do
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

      ! get local indices of elements distributed to/from other processors
      nn_to = 0
      nn_from = 0

      ! loop over all vmu indices, find corresponding y indices
      do ikxyz = kxyz_lo%llim_world, kxyz_lo%ulim_world
         do imu = 1, nmu
            do iv = 1, nvpa
               ! obtain corresponding y indices
               call kxyzidx2vmuidx(iv, imu, ikxyz, kxyz_lo, vmu_lo, iy, ikx, iz, it, ivmu)
               ! if vmu index local, set:
               ! ip = corresponding y processor
               ! from_list%first-third arrays = iv,imu,ikxyz  (ie vmu indices)
               ! later will send from_list to proc ip
               if (idx_local(kxyz_lo, ikxyz)) then
                  ip = proc_id(vmu_lo, ivmu)
                  n = nn_from(ip) + 1
                  nn_from(ip) = n
                  from_list(ip)%first(n) = iv
                  from_list(ip)%second(n) = imu
                  from_list(ip)%third(n) = ikxyz
               end if
               ! if y index local, set ip to corresponding vmu processor
               ! set to_list%first,second arrays = iy,iy  (ie y indices)
               ! will receive to_list from ip
               if (idx_local(vmu_lo, ivmu)) then
                  ip = proc_id(kxyz_lo, ikxyz)
                  n = nn_to(ip) + 1
                  nn_to(ip) = n
                  to_list(ip)%first(n) = iy
                  to_list(ip)%second(n) = ikx
                  to_list(ip)%third(n) = iz
                  to_list(ip)%fourth(n) = it
                  to_list(ip)%fifth(n) = ivmu
               end if
            end do
         end do
      end do

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

      call init_redist(kxyz2vmu, 'c', to_low, to_high, to_list, &
                       from_low, from_high, from_list)

      call delete_list(to_list)
      call delete_list(from_list)

   end subroutine init_kxyz_to_vmu_redistribute

   subroutine init_xyz_to_vmu_redistribute

      use mp, only: nproc
      use stella_layouts, only: xyz_lo, vmu_lo
      use stella_layouts, only: xyzidx2vmuidx
      use stella_layouts, only: idx_local, proc_id
      use redistribute, only: index_list_type, init_redist
      use redistribute, only: delete_list, set_redist_character_type
      use grids_velocity, only: nvpa, nmu
      use grids_z, only: nzgrid

      implicit none

      type(index_list_type), dimension(0:nproc - 1) :: to_list, from_list
      integer, dimension(0:nproc - 1) :: nn_to, nn_from
      integer, dimension(3) :: from_low, from_high
      integer, dimension(5) :: to_high, to_low
      integer :: ixyz, ivmu
      integer :: iv, imu, iy, ix, iz, it
      integer :: ip, n
      logical :: initialized = .false.

      if (initialized) return
      initialized = .true.

      ! count number of elements to be redistributed to/from each processor
      nn_to = 0
      nn_from = 0
      do ixyz = xyz_lo%llim_world, xyz_lo%ulim_world
         do imu = 1, nmu
            do iv = 1, nvpa
               call xyzidx2vmuidx(iv, imu, ixyz, xyz_lo, vmu_lo, iy, ix, iz, it, ivmu)
               if (idx_local(xyz_lo, ixyz)) &
                  nn_from(proc_id(vmu_lo, ivmu)) = nn_from(proc_id(vmu_lo, ivmu)) + 1
               if (idx_local(vmu_lo, ivmu)) &
                  nn_to(proc_id(xyz_lo, ixyz)) = nn_to(proc_id(xyz_lo, ixyz)) + 1
            end do
         end do
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

      ! get local indices of elements distributed to/from other processors
      nn_to = 0
      nn_from = 0

      ! loop over all vmu indices, find corresponding y indices
      do ixyz = xyz_lo%llim_world, xyz_lo%ulim_world
         do imu = 1, nmu
            do iv = 1, nvpa
               ! obtain corresponding y indices
               call xyzidx2vmuidx(iv, imu, ixyz, xyz_lo, vmu_lo, iy, ix, iz, it, ivmu)
               ! if vmu index local, set:
               ! ip = corresponding y processor
               ! from_list%first-third arrays = iv,imu,ixyz  (ie vmu indices)
               ! later will send from_list to proc ip
               if (idx_local(xyz_lo, ixyz)) then
                  ip = proc_id(vmu_lo, ivmu)
                  n = nn_from(ip) + 1
                  nn_from(ip) = n
                  from_list(ip)%first(n) = iv
                  from_list(ip)%second(n) = imu
                  from_list(ip)%third(n) = ixyz
               end if
               ! if y index local, set ip to corresponding vmu processor
               ! set to_list%first,second arrays = iy,iy  (ie y indices)
               ! will receive to_list from ip
               if (idx_local(vmu_lo, ivmu)) then
                  ip = proc_id(xyz_lo, ixyz)
                  n = nn_to(ip) + 1
                  nn_to(ip) = n
                  to_list(ip)%first(n) = iy
                  to_list(ip)%second(n) = ix
                  to_list(ip)%third(n) = iz
                  to_list(ip)%fourth(n) = it
                  to_list(ip)%fifth(n) = ivmu
               end if
            end do
         end do
      end do

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

      call init_redist(xyz2vmu, 'r', to_low, to_high, to_list, &
                       from_low, from_high, from_list)

      call delete_list(to_list)
      call delete_list(from_list)

   end subroutine init_xyz_to_vmu_redistribute

   subroutine init_kymus_to_vmus_redistribute

      use mp, only: nproc
      use stella_layouts, only: kymus_lo, vmu_lo
      use stella_layouts, only: kymusidx2vmuidx
      use stella_layouts, only: idx_local, proc_id
      use redistribute, only: index_list_type, init_redist
      use redistribute, only: delete_list, set_redist_character_type
      use grids_velocity, only: nvpa
      use parameters_kxky_grid, only: nakx
      use grids_z, only: nzgrid, ntubes

      implicit none

      type(index_list_type), dimension(0:nproc - 1) :: to_list, from_list
      integer, dimension(0:nproc - 1) :: nn_to, nn_from
      integer, dimension(5) :: from_low, from_high
      integer, dimension(5) :: to_high, to_low
      integer :: ikymus, ivmu
      integer :: iv, iky, ikx, iz, it
      integer :: ip, n
      logical :: initialized = .false.

      if (initialized) return
      initialized = .true.

      ! count number of elements to be redistributed to/from each processor
      nn_to = 0
      nn_from = 0
      do ikymus = kymus_lo%llim_world, kymus_lo%ulim_world
         do iv = 1, nvpa
            do it = 1, ntubes
               do iz = -nzgrid, nzgrid
                  do ikx = 1, nakx
                     call kymusidx2vmuidx(iv, ikymus, kymus_lo, vmu_lo, iky, ivmu)
                     if (idx_local(kymus_lo, ikymus)) &
                          nn_from(proc_id(vmu_lo, ivmu)) = nn_from(proc_id(vmu_lo, ivmu)) + 1
                     if (idx_local(vmu_lo, ivmu)) &
                          nn_to(proc_id(kymus_lo, ikymus)) = nn_to(proc_id(kymus_lo, ikymus)) + 1
                  end do
               end do
            end do
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

      ! get local indices of elements distributed to/from other processors
      nn_to = 0
      nn_from = 0

      ! loop over all indices in the kymus layout and find the corresponding indices for iky and ivmu
      do ikymus = kymus_lo%llim_world, kymus_lo%ulim_world
         do iv = 1, nvpa
            do it = 1, ntubes
               do iz = -nzgrid, nzgrid
                  do ikx = 1, nakx
                     ! obtain corresponding ky indices
                     call kymusidx2vmuidx(iv, ikymus, kymus_lo, vmu_lo, iky, ivmu)
                     ! if kymus index local, set:
                     ! ip = corresponding y processor
                     ! from_list%first-third arrays = ikx,iz,it,iv,ikymus
                     ! later will send from_list to proc ip
                     if (idx_local(kymus_lo, ikymus)) then
                        ip = proc_id(vmu_lo, ivmu)
                        n = nn_from(ip) + 1
                        nn_from(ip) = n
                        from_list(ip)%first(n) = ikx
                        from_list(ip)%second(n) = iz
                        from_list(ip)%third(n) = it
                        from_list(ip)%fourth(n) = iv
                        from_list(ip)%fifth(n) = ikymus
                     end if
                     if (idx_local(vmu_lo, ivmu)) then
                        ip = proc_id(kymus_lo, ikymus)
                        n = nn_to(ip) + 1
                        nn_to(ip) = n
                        to_list(ip)%first(n) = iky
                        to_list(ip)%second(n) = ikx
                        to_list(ip)%third(n) = iz
                        to_list(ip)%fourth(n) = it
                        to_list(ip)%fifth(n) = ivmu
                     end if
                  end do
               end do
            end do
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

      call init_redist(kymus2vmus, 'c', to_low, to_high, to_list, &
                       from_low, from_high, from_list)

      call delete_list(to_list)
      call delete_list(from_list)

   end subroutine init_kymus_to_vmus_redistribute

   subroutine test_kymus_to_vmus_redistribute

      use redistribute, only: scatter, gather
      use store_arrays_distribution_fn, only: g_kymus, gnew
      use grids_velocity, only: nvpa
      use mp, only: proc0, send, receive
      use grids_z, only: nzgrid, ntubes
      use parameters_kxky_grid, only: nakx, naky
      use stella_layouts, only: vmu_lo, kymus_lo
      use stella_layouts, only: iv_idx, imu_idx, is_idx, iky_idx, idx_local, proc_id

      implicit none

      complex, dimension (:, :, :, :), allocatable :: gtmp
      integer :: ikymus, ivmu, iv, imu, is, it, iz, ikx, iky, izmod

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
      deallocate (gtmp)

      call scatter(kymus2vmus, gnew, g_kymus)
      
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
      deallocate (gtmp)

      call gather(kymus2vmus, g_kymus, gnew)
      
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
      deallocate (gtmp)
     
   end subroutine test_kymus_to_vmus_redistribute
   
   subroutine finish_redistribute

      implicit none

      redistribute_initialized = .false.

   end subroutine finish_redistribute

end module calculations_redistribute
