module stella_layouts

  use common_types, only: vmu_layout_type
  use common_types, only: kxkyz_layout_type, kxyz_layout_type, xyz_layout_type
  
  implicit none

  private

  public :: xyzs_layout, vms_layout
  public :: finish_layouts
  
  public :: init_stella_layouts, init_dist_fn_layouts
  public :: kxkyz_lo, kxyz_lo, xyz_lo, vmu_lo
  
  public :: kxkyzidx2vmuidx, kxyzidx2vmuidx, xyzidx2vmuidx
  
  public :: iz_idx, iky_idx, ikx_idx, iv_idx, imu_idx, is_idx, iy_idx
  public :: idx, proc_id, idx_local
  
  character (len=4) :: xyzs_layout
  character (len=3) :: vms_layout
  logical :: exist
  
  type (kxkyz_layout_type) :: kxkyz_lo
  type (kxyz_layout_type) :: kxyz_lo
  type (xyz_layout_type) :: xyz_lo
  type (vmu_layout_type) :: vmu_lo
  
  interface iz_idx
     module procedure iz_idx_kxkyz
     module procedure iz_idx_kxyz
     module procedure iz_idx_xyz
  end interface

  interface iv_idx
     module procedure iv_idx_vmu
  end interface

  interface iky_idx
     module procedure iky_idx_kxkyz
  end interface

  interface iy_idx
     module procedure iy_idx_kxyz
     module procedure iy_idx_xyz
  end interface

  interface ikx_idx
     module procedure ikx_idx_kxkyz
     module procedure ikx_idx_kxyz
  end interface

  interface ix_idx
     module procedure ix_idx_xyz
  end interface

  interface imu_idx
     module procedure imu_idx_vmu
  end interface

  interface is_idx
     module procedure is_idx_kxkyz
     module procedure is_idx_kxyz
     module procedure is_idx_xyz
     module procedure is_idx_vmu
  end interface

  interface proc_id
     module procedure proc_id_kxkyz
     module procedure proc_id_kxyz
     module procedure proc_id_xyz
     module procedure proc_id_vmu
  end interface

  interface idx
     module procedure idx_kxkyz
     module procedure idx_kxyz
     module procedure idx_xyz
     module procedure idx_vmu
  end interface

  interface idx_local
     module procedure idx_local_kxkyz, ig_local_kxkyz
     module procedure idx_local_kxyz, ig_local_kxyz
     module procedure idx_local_xyz, ig_local_xyz
     module procedure idx_local_vmu, ig_local_vmu
  end interface

contains

  subroutine init_stella_layouts
    
    use mp, only: proc0
    implicit none
    logical, save :: initialized = .false.

    if (initialized) return
    initialized = .true.
    
    if (proc0) call read_parameters
    call broadcast_results

  end subroutine init_stella_layouts

  subroutine read_parameters

    use mp, only: mp_abort
    use file_utils, only: input_unit, error_unit, input_unit_exist, error_unit

    implicit none

    integer :: in_file

    namelist /layouts_knobs/ xyzs_layout, vms_layout

    xyzs_layout = 'xyzs'
    vms_layout = 'vms'

    in_file=input_unit_exist("layouts_knobs", exist)
    if (exist) read (unit=input_unit("layouts_knobs"), nml=layouts_knobs)

    if (xyzs_layout.ne.'xyzs' .and. &
         xyzs_layout.ne.'xzys' .and. &
         xyzs_layout.ne.'yxzs' .and. &
         xyzs_layout.ne.'yzxs' .and. &
         xyzs_layout.ne.'zxys' .and. &
         xyzs_layout.ne.'zyxs') then
       call mp_abort ('stella_layouts: read_parameters finds illegal xyzs_layout. aborting')
    endif
    if (vms_layout.ne.'vms' .and. &
         vms_layout.ne.'mvs') then
       call mp_abort ('stella_layouts: read_parameters finds illegal vms_layout. aborting')
    end if

  end subroutine read_parameters
    
  subroutine broadcast_results
    use mp, only: broadcast
    implicit none
    call broadcast (xyzs_layout)
    call broadcast (vms_layout)
  end subroutine broadcast_results

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Distribution function layouts
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine init_dist_fn_layouts (nzgrid, naky, nakx, nvgrid, nmu, nspec, ny, nx)

    implicit none

    integer, intent (in) :: nzgrid, naky, nakx, nvgrid, nmu, nspec, ny, nx

    call init_kxkyz_layout (nzgrid, naky, nakx, nvgrid, nmu, nspec)
    call init_kxyz_layout (nzgrid, naky, nakx, nvgrid, nmu, nspec, ny)
    call init_xyz_layout (nzgrid, naky, nakx, nvgrid, nmu, nspec, ny, nx)
    call init_vmu_layout (nzgrid, naky, nakx, nvgrid, nmu, nspec, ny, nx)

  end subroutine init_dist_fn_layouts

  subroutine init_kxkyz_layout &
       (nzgrid, naky, nakx, nvgrid, nmu, nspec)
    
    use mp, only: iproc, nproc

    implicit none

    integer, intent (in) :: nzgrid, naky, nakx, nvgrid, nmu, nspec
    logical, save :: initialized = .false.

    if (initialized) return
    initialized = .true.
   
    kxkyz_lo%iproc = iproc
    kxkyz_lo%nzgrid = nzgrid
    kxkyz_lo%nzed = 2*nzgrid+1
    kxkyz_lo%naky = naky
    kxkyz_lo%nakx = nakx
    kxkyz_lo%nvgrid = nvgrid
    kxkyz_lo%nvpa = 2*nvgrid+1
    kxkyz_lo%nmu = nmu
    kxkyz_lo%nspec = nspec
    kxkyz_lo%llim_world = 0
    kxkyz_lo%ulim_world = naky*nakx*kxkyz_lo%nzed*nspec - 1
    kxkyz_lo%blocksize = kxkyz_lo%ulim_world/nproc + 1
    kxkyz_lo%llim_proc = kxkyz_lo%blocksize*iproc
    kxkyz_lo%ulim_proc = min(kxkyz_lo%ulim_world, kxkyz_lo%llim_proc + kxkyz_lo%blocksize - 1)
    kxkyz_lo%ulim_alloc = max(kxkyz_lo%llim_proc, kxkyz_lo%ulim_proc)

  end subroutine init_kxkyz_layout

  elemental function is_idx_kxkyz (lo, i)
    implicit none
    integer :: is_idx_kxkyz
    type (kxkyz_layout_type), intent (in) :: lo
    integer, intent (in) :: i
    is_idx_kxkyz = 1 + mod((i - lo%llim_world)/lo%nakx/lo%naky/lo%nzed, lo%nspec)
  end function is_idx_kxkyz

  elemental function ikx_idx_kxkyz (lo, i)

    implicit none

    integer :: ikx_idx_kxkyz
    type (kxkyz_layout_type), intent (in) :: lo
    integer, intent (in) :: i

    select case (xyzs_layout)
    case ('xyzs')
       ikx_idx_kxkyz = 1 + mod((i - lo%llim_world), lo%nakx)
    case ('xzys')
       ikx_idx_kxkyz = 1 + mod((i - lo%llim_world), lo%nakx)
    case ('yxzs')
       ikx_idx_kxkyz = 1 + mod((i - lo%llim_world)/lo%naky, lo%nakx)
    case ('zxys')
       ikx_idx_kxkyz = 1 + mod((i - lo%llim_world)/lo%nzed, lo%nakx)
    case ('zyxs')
       ikx_idx_kxkyz = 1 + mod((i - lo%llim_world)/lo%nzed/lo%naky, lo%nakx)
    case ('yzxs')
       ikx_idx_kxkyz = 1 + mod((i - lo%llim_world)/lo%naky/lo%nzed, lo%nakx)
    end select

  end function ikx_idx_kxkyz

  elemental function iky_idx_kxkyz (lo, i)

    implicit none
    integer :: iky_idx_kxkyz
    type (kxkyz_layout_type), intent (in) :: lo
    integer, intent (in) :: i

    select case (xyzs_layout)
    case ('yxzs')
       iky_idx_kxkyz = 1 + mod(i - lo%llim_world, lo%naky)
    case ('yzxs')
       iky_idx_kxkyz = 1 + mod(i - lo%llim_world, lo%naky)
    case ('xyzs')
       iky_idx_kxkyz = 1 + mod((i - lo%llim_world)/lo%nakx, lo%naky)
    case ('zyxs')
       iky_idx_kxkyz = 1 + mod((i - lo%llim_world)/lo%nzed, lo%naky)
    case ('zxys')
       iky_idx_kxkyz = 1 + mod((i - lo%llim_world)/lo%nzed/lo%nakx, lo%naky)
    case ('xzys')
       iky_idx_kxkyz = 1 + mod((i - lo%llim_world)/lo%nakx/lo%nzed, lo%naky)
    end select

  end function iky_idx_kxkyz

  elemental function iz_idx_kxkyz (lo, i)

    implicit none
    integer :: iz_idx_kxkyz
    type (kxkyz_layout_type), intent (in) :: lo
    integer, intent (in) :: i

    select case (xyzs_layout)
    case ('zyxs')
       iz_idx_kxkyz = -lo%nzgrid + mod((i - lo%llim_world), lo%nzed)
    case ('zxys')
       iz_idx_kxkyz = -lo%nzgrid + mod((i - lo%llim_world), lo%nzed)
    case ('xzys')
       iz_idx_kxkyz = -lo%nzgrid + mod((i - lo%llim_world)/lo%nakx, lo%nzed)
    case ('yzxs')
       iz_idx_kxkyz = -lo%nzgrid + mod((i - lo%llim_world)/lo%naky, lo%nzed)
    case ('yxzs')
       iz_idx_kxkyz = -lo%nzgrid + mod((i - lo%llim_world)/lo%naky/lo%nakx, lo%nzed)
    case ('xyzs')
       iz_idx_kxkyz = -lo%nzgrid + mod((i - lo%llim_world)/lo%nakx/lo%naky, lo%nzed)
    end select

  end function iz_idx_kxkyz

  elemental function proc_id_kxkyz (lo, i)
    implicit none
    integer :: proc_id_kxkyz
    type (kxkyz_layout_type), intent (in) :: lo
    integer, intent (in) :: i

    proc_id_kxkyz = i/lo%blocksize

  end function proc_id_kxkyz

  elemental function idx_kxkyz (lo, iky, ikx, ig, is)

    implicit none

    integer :: idx_kxkyz
    type (kxkyz_layout_type), intent (in) :: lo
    integer, intent (in) :: iky, ikx, ig, is

    select case (xyzs_layout)
    case ('xyzs')
       idx_kxkyz = ikx-1 + lo%nakx*(iky-1 + lo%naky*(ig+lo%nzgrid + lo%nzed*(is-1)))
    case ('xzys')
       idx_kxkyz = ikx-1 + lo%nakx*(ig+lo%nzgrid + lo%nzed*(iky-1 + lo%naky*(is-1)))
    case ('yxzs')
       idx_kxkyz = iky-1 + lo%naky*(ikx-1 + lo%nakx*(ig+lo%nzgrid + lo%nzed*(is-1)))
    case ('yzxs')
       idx_kxkyz = iky-1 + lo%naky*(ig+lo%nzgrid + lo%nzed*(ikx-1 + lo%nakx*(is-1)))
    case ('zyxs')
       idx_kxkyz = ig+lo%nzgrid + lo%nzed*(iky-1 + lo%naky*(ikx-1))
    case ('zxys')
       idx_kxkyz = ig+lo%nzgrid + lo%nzed*(ikx-1 + lo%nakx*(iky-1))
    end select

  end function idx_kxkyz

  elemental function idx_local_kxkyz (lo, iky, ikx, ig, is)

    implicit none
    logical :: idx_local_kxkyz
    type (kxkyz_layout_type), intent (in) :: lo
    integer, intent (in) :: iky, ikx, ig, is

    idx_local_kxkyz = idx_local(lo, idx(lo, iky, ikx, ig, is))
  end function idx_local_kxkyz

  elemental function ig_local_kxkyz (lo, ig)
    implicit none
    logical :: ig_local_kxkyz
    type (kxkyz_layout_type), intent (in) :: lo
    integer, intent (in) :: ig

    ig_local_kxkyz = lo%iproc == proc_id(lo, ig)
  end function ig_local_kxkyz

  subroutine init_kxyz_layout &
       (nzgrid, naky, nakx, nvgrid, nmu, nspec, ny)
    
    use mp, only: iproc, nproc

    implicit none

    integer, intent (in) :: nzgrid, ny, naky, nakx, nvgrid, nmu, nspec
    logical, save :: initialized = .false.

    if (initialized) return
    initialized = .true.
   
    kxyz_lo%iproc = iproc
    kxyz_lo%nzgrid = nzgrid
    kxyz_lo%nzed = 2*nzgrid+1
    kxyz_lo%ny = ny
    kxyz_lo%naky = naky
    kxyz_lo%nakx = nakx
    kxyz_lo%nvgrid = nvgrid
    kxyz_lo%nvpa = 2*nvgrid+1
    kxyz_lo%nmu = nmu
    kxyz_lo%nspec = nspec
    kxyz_lo%llim_world = 0
    kxyz_lo%ulim_world = ny*nakx*kxyz_lo%nzed*nspec - 1
    kxyz_lo%blocksize = kxyz_lo%ulim_world/nproc + 1
    kxyz_lo%llim_proc = kxyz_lo%blocksize*iproc
    kxyz_lo%ulim_proc = min(kxyz_lo%ulim_world, kxyz_lo%llim_proc + kxyz_lo%blocksize - 1)
    kxyz_lo%ulim_alloc = max(kxyz_lo%llim_proc, kxyz_lo%ulim_proc)

  end subroutine init_kxyz_layout

  elemental function is_idx_kxyz (lo, i)
    implicit none
    integer :: is_idx_kxyz
    type (kxyz_layout_type), intent (in) :: lo
    integer, intent (in) :: i
    is_idx_kxyz = 1 + mod((i - lo%llim_world)/lo%nakx/lo%ny/lo%nzed, lo%nspec)
  end function is_idx_kxyz

  elemental function ikx_idx_kxyz (lo, i)

    implicit none

    integer :: ikx_idx_kxyz
    type (kxyz_layout_type), intent (in) :: lo
    integer, intent (in) :: i

    select case (xyzs_layout)
    case ('xyzs')
       ikx_idx_kxyz = 1 + mod((i - lo%llim_world), lo%nakx)
    case ('xzys')
       ikx_idx_kxyz = 1 + mod((i - lo%llim_world), lo%nakx)
    case ('yxzs')
       ikx_idx_kxyz = 1 + mod((i - lo%llim_world)/lo%ny, lo%nakx)
    case ('zxys')
       ikx_idx_kxyz = 1 + mod((i - lo%llim_world)/lo%nzed, lo%nakx)
    case ('zyxs')
       ikx_idx_kxyz = 1 + mod((i - lo%llim_world)/lo%nzed/lo%ny, lo%nakx)
    case ('yzxs')
       ikx_idx_kxyz = 1 + mod((i - lo%llim_world)/lo%ny/lo%nzed, lo%nakx)
    end select

  end function ikx_idx_kxyz

  elemental function iy_idx_kxyz (lo, i)

    implicit none
    integer :: iy_idx_kxyz
    type (kxyz_layout_type), intent (in) :: lo
    integer, intent (in) :: i

    select case (xyzs_layout)
    case ('yxzs')
       iy_idx_kxyz = 1 + mod(i - lo%llim_world, lo%ny)
    case ('yzxs')
       iy_idx_kxyz = 1 + mod(i - lo%llim_world, lo%ny)
    case ('xyzs')
       iy_idx_kxyz = 1 + mod((i - lo%llim_world)/lo%nakx, lo%ny)
    case ('zyxs')
       iy_idx_kxyz = 1 + mod((i - lo%llim_world)/lo%nzed, lo%ny)
    case ('zxys')
       iy_idx_kxyz = 1 + mod((i - lo%llim_world)/lo%nzed/lo%nakx, lo%ny)
    case ('xzys')
       iy_idx_kxyz = 1 + mod((i - lo%llim_world)/lo%nakx/lo%nzed, lo%ny)
    end select

  end function iy_idx_kxyz

  elemental function iz_idx_kxyz (lo, i)

    implicit none
    integer :: iz_idx_kxyz
    type (kxyz_layout_type), intent (in) :: lo
    integer, intent (in) :: i

    select case (xyzs_layout)
    case ('zyxs')
       iz_idx_kxyz = -lo%nzgrid + mod((i - lo%llim_world), lo%nzed)
    case ('zxys')
       iz_idx_kxyz = -lo%nzgrid + mod((i - lo%llim_world), lo%nzed)
    case ('yzxs')
       iz_idx_kxyz = -lo%nzgrid + mod((i - lo%llim_world)/lo%ny, lo%nzed)
    case ('xzys')
       iz_idx_kxyz = -lo%nzgrid + mod((i - lo%llim_world)/lo%nakx, lo%nzed)
    case ('yxzs')
       iz_idx_kxyz = -lo%nzgrid + mod((i - lo%llim_world)/lo%ny/lo%nakx, lo%nzed)
    case ('xyzs')
       iz_idx_kxyz = -lo%nzgrid + mod((i - lo%llim_world)/lo%nakx/lo%ny, lo%nzed)
    end select

  end function iz_idx_kxyz

  elemental function proc_id_kxyz (lo, i)
    implicit none
    integer :: proc_id_kxyz
    type (kxyz_layout_type), intent (in) :: lo
    integer, intent (in) :: i

    proc_id_kxyz = i/lo%blocksize

  end function proc_id_kxyz

  elemental function idx_kxyz (lo, iy, ikx, ig, is)

    implicit none

    integer :: idx_kxyz
    type (kxyz_layout_type), intent (in) :: lo
    integer, intent (in) :: iy, ikx, ig, is

    select case (xyzs_layout)
    case ('xyzs')
       idx_kxyz = ikx-1 + lo%nakx*(iy-1 + lo%ny*(ig+lo%nzgrid + lo%nzed*(is-1)))
    case ('xzys')
       idx_kxyz = ikx-1 + lo%nakx*(ig+lo%nzgrid + lo%nzed*(iy-1 + lo%ny*(is-1)))
    case ('yxzs')
       idx_kxyz = iy-1 + lo%ny*(ikx-1 + lo%nakx*(ig+lo%nzgrid + lo%nzed*(is-1)))
    case ('yzxs')
       idx_kxyz = iy-1 + lo%ny*(ig+lo%nzgrid + lo%nzed*(ikx-1 + lo%nakx*(is-1)))
    case ('zyxs')
       idx_kxyz = ig+lo%nzgrid + lo%nzed*(iy-1 + lo%ny*(ikx-1))
    case ('zxys')
       idx_kxyz = ig+lo%nzgrid + lo%nzed*(ikx-1 + lo%nakx*(iy-1))
    end select

  end function idx_kxyz

  elemental function idx_local_kxyz (lo, iy, ikx, ig, is)

    implicit none
    logical :: idx_local_kxyz
    type (kxyz_layout_type), intent (in) :: lo
    integer, intent (in) :: iy, ikx, ig, is

    idx_local_kxyz = idx_local(lo, idx(lo, iy, ikx, ig, is))
  end function idx_local_kxyz

  elemental function ig_local_kxyz (lo, ig)
    implicit none
    logical :: ig_local_kxyz
    type (kxyz_layout_type), intent (in) :: lo
    integer, intent (in) :: ig

    ig_local_kxyz = lo%iproc == proc_id(lo, ig)
  end function ig_local_kxyz

  subroutine init_xyz_layout &
       (nzgrid, naky, nakx, nvgrid, nmu, nspec, ny, nx)
    
    use mp, only: iproc, nproc

    implicit none

    integer, intent (in) :: nzgrid, ny, nx, naky, nakx, nvgrid, nmu, nspec
    logical, save :: initialized = .false.

    if (initialized) return
    initialized = .true.
   
    xyz_lo%iproc = iproc
    xyz_lo%nzgrid = nzgrid
    xyz_lo%nzed = 2*nzgrid+1
    xyz_lo%ny = ny
    xyz_lo%nx = nx
    xyz_lo%naky = naky
    xyz_lo%nakx = nakx
    xyz_lo%nvgrid = nvgrid
    xyz_lo%nvpa = 2*nvgrid+1
    xyz_lo%nmu = nmu
    xyz_lo%nspec = nspec
    xyz_lo%llim_world = 0
    xyz_lo%ulim_world = ny*nx*xyz_lo%nzed*nspec - 1
    xyz_lo%blocksize = xyz_lo%ulim_world/nproc + 1
    xyz_lo%llim_proc = xyz_lo%blocksize*iproc
    xyz_lo%ulim_proc = min(xyz_lo%ulim_world, xyz_lo%llim_proc + xyz_lo%blocksize - 1)
    xyz_lo%ulim_alloc = max(xyz_lo%llim_proc, xyz_lo%ulim_proc)

  end subroutine init_xyz_layout

  elemental function is_idx_xyz (lo, i)
    implicit none
    integer :: is_idx_xyz
    type (xyz_layout_type), intent (in) :: lo
    integer, intent (in) :: i
    is_idx_xyz = 1 + mod((i - lo%llim_world)/lo%nx/lo%ny/lo%nzed, lo%nspec)
  end function is_idx_xyz

  elemental function ix_idx_xyz (lo, i)

    implicit none

    integer :: ix_idx_xyz
    type (xyz_layout_type), intent (in) :: lo
    integer, intent (in) :: i

    select case (xyzs_layout)
    case ('xyzs')
       ix_idx_xyz = 1 + mod((i - lo%llim_world), lo%nx)
    case ('xzys')
       ix_idx_xyz = 1 + mod((i - lo%llim_world), lo%nx)
    case ('yxzs')
       ix_idx_xyz = 1 + mod((i - lo%llim_world)/lo%ny, lo%nx)
    case ('zxys')
       ix_idx_xyz = 1 + mod((i - lo%llim_world)/lo%nzed, lo%nx)
    case ('zyxs')
       ix_idx_xyz = 1 + mod((i - lo%llim_world)/lo%nzed/lo%ny, lo%nx)
    case ('yzxs')
       ix_idx_xyz = 1 + mod((i - lo%llim_world)/lo%ny/lo%nzed, lo%nx)
    end select

  end function ix_idx_xyz

  elemental function iy_idx_xyz (lo, i)

    implicit none
    integer :: iy_idx_xyz
    type (xyz_layout_type), intent (in) :: lo
    integer, intent (in) :: i

    select case (xyzs_layout)
    case ('yxzs')
       iy_idx_xyz = 1 + mod(i - lo%llim_world, lo%ny)
    case ('yzxs')
       iy_idx_xyz = 1 + mod(i - lo%llim_world, lo%ny)
    case ('xyzs')
       iy_idx_xyz = 1 + mod((i - lo%llim_world)/lo%nx, lo%ny)
    case ('zyxs')
       iy_idx_xyz = 1 + mod((i - lo%llim_world)/lo%nzed, lo%ny)
    case ('zxys')
       iy_idx_xyz = 1 + mod((i - lo%llim_world)/lo%nzed/lo%nx, lo%ny)
    case ('xzys')
       iy_idx_xyz = 1 + mod((i - lo%llim_world)/lo%nx/lo%nzed, lo%ny)
    end select

  end function iy_idx_xyz

  elemental function iz_idx_xyz (lo, i)

    implicit none
    integer :: iz_idx_xyz
    type (xyz_layout_type), intent (in) :: lo
    integer, intent (in) :: i

    select case (xyzs_layout)
    case ('zyxs')
       iz_idx_xyz = -lo%nzgrid + mod((i - lo%llim_world), lo%nzed)
    case ('zxys')
       iz_idx_xyz = -lo%nzgrid + mod((i - lo%llim_world), lo%nzed)
    case ('yzxs')
       iz_idx_xyz = -lo%nzgrid + mod((i - lo%llim_world)/lo%ny, lo%nzed)
    case ('xzys')
       iz_idx_xyz = -lo%nzgrid + mod((i - lo%llim_world)/lo%nx, lo%nzed)
    case ('yxzs')
       iz_idx_xyz = -lo%nzgrid + mod((i - lo%llim_world)/lo%ny/lo%nx, lo%nzed)
    case ('xyzs')
       iz_idx_xyz = -lo%nzgrid + mod((i - lo%llim_world)/lo%nx/lo%ny, lo%nzed)
    end select

  end function iz_idx_xyz

  elemental function proc_id_xyz (lo, i)
    implicit none
    integer :: proc_id_xyz
    type (xyz_layout_type), intent (in) :: lo
    integer, intent (in) :: i

    proc_id_xyz = i/lo%blocksize

  end function proc_id_xyz

  elemental function idx_xyz (lo, iy, ix, ig, is)

    implicit none

    integer :: idx_xyz
    type (xyz_layout_type), intent (in) :: lo
    integer, intent (in) :: iy, ix, ig, is

    select case (xyzs_layout)
    case ('xyzs')
       idx_xyz = ix-1 + lo%nx*(iy-1 + lo%ny*(ig+lo%nzgrid + lo%nzed*(is-1)))
    case ('xzys')
       idx_xyz = ix-1 + lo%nx*(ig+lo%nzgrid + lo%nzed*(iy-1 + lo%ny*(is-1)))
    case ('yxzs')
       idx_xyz = iy-1 + lo%ny*(ix-1 + lo%nx*(ig+lo%nzgrid + lo%nzed*(is-1)))
    case ('yzxs')
       idx_xyz = iy-1 + lo%ny*(ig+lo%nzgrid + lo%nzed*(ix-1 + lo%nx*(is-1)))
    case ('zyxs')
       idx_xyz = ig+lo%nzgrid + lo%nzed*(iy-1 + lo%ny*(ix-1))
    case ('zxys')
       idx_xyz = ig+lo%nzgrid + lo%nzed*(ix-1 + lo%nx*(iy-1))
    end select

  end function idx_xyz

  elemental function idx_local_xyz (lo, iy, ix, ig, is)

    implicit none
    logical :: idx_local_xyz
    type (xyz_layout_type), intent (in) :: lo
    integer, intent (in) :: iy, ix, ig, is

    idx_local_xyz = idx_local(lo, idx(lo, iy, ix, ig, is))
  end function idx_local_xyz

  elemental function ig_local_xyz (lo, ig)
    implicit none
    logical :: ig_local_xyz
    type (xyz_layout_type), intent (in) :: lo
    integer, intent (in) :: ig

    ig_local_xyz = lo%iproc == proc_id(lo, ig)
  end function ig_local_xyz

  subroutine init_vmu_layout &
       (nzgrid, naky, nakx, nvgrid, nmu, nspec, ny, nx)
    
    use mp, only: iproc, nproc

    implicit none

    integer, intent (in) :: nzgrid, naky, nakx, nvgrid, nmu, nspec, ny, nx
    logical, save :: initialized = .false.

    if (initialized) return
    initialized = .true.
   
    vmu_lo%xyz = .true.
    vmu_lo%iproc = iproc
    vmu_lo%nzed = 2*nzgrid+1
    vmu_lo%nzgrid = nzgrid
    vmu_lo%ny = ny
    vmu_lo%naky = naky
    vmu_lo%nx = nx
    vmu_lo%nakx = nakx
    vmu_lo%nvgrid = nvgrid
    vmu_lo%nvpa = 2*nvgrid+1
    vmu_lo%nmu = nmu
    vmu_lo%nspec = nspec
    vmu_lo%llim_world = 0
    vmu_lo%ulim_world = vmu_lo%nvpa*nmu*nspec - 1
      
    vmu_lo%blocksize = vmu_lo%ulim_world/nproc + 1
    vmu_lo%llim_proc = vmu_lo%blocksize*iproc
    vmu_lo%ulim_proc = min(vmu_lo%ulim_world, vmu_lo%llim_proc + vmu_lo%blocksize - 1)
    vmu_lo%ulim_alloc = max(vmu_lo%llim_proc, vmu_lo%ulim_proc)

  end subroutine init_vmu_layout

  elemental function is_idx_vmu (lo, i)

    implicit none
    integer :: is_idx_vmu
    type (vmu_layout_type), intent (in) :: lo
    integer, intent (in) :: i

    ! the order of the division does not matter, so no need for branching
    is_idx_vmu = 1 + mod((i - lo%llim_world)/lo%nvpa/lo%nmu, lo%nspec)

  end function is_idx_vmu

  elemental function imu_idx_vmu (lo, i)

    implicit none

    integer :: imu_idx_vmu
    type (vmu_layout_type), intent (in) :: lo
    integer, intent (in) :: i

    select case (vms_layout)
    case ('vms')
       imu_idx_vmu = 1 + mod((i - lo%llim_world)/lo%nvpa, lo%nmu)
    case ('mvs')
       imu_idx_vmu = 1 + mod((i - lo%llim_world), lo%nmu)
    end select

  end function imu_idx_vmu

  elemental function iv_idx_vmu (lo, i)

    implicit none
    integer :: iv_idx_vmu
    type (vmu_layout_type), intent (in) :: lo
    integer, intent (in) :: i

    select case (vms_layout)
    case ('vms')
       iv_idx_vmu = -lo%nvgrid + mod((i - lo%llim_world), lo%nvpa)
    case ('mvs')
       iv_idx_vmu = -lo%nvgrid + mod((i - lo%llim_world)/lo%nmu, lo%nvpa)
    end select

  end function iv_idx_vmu

  elemental function proc_id_vmu (lo, i)
    implicit none
    integer :: proc_id_vmu
    type (vmu_layout_type), intent (in) :: lo
    integer, intent (in) :: i

    proc_id_vmu = i/lo%blocksize

  end function proc_id_vmu

  elemental function idx_vmu (lo, iv, imu, is)

    implicit none

    integer :: idx_vmu
    type (vmu_layout_type), intent (in) :: lo
    integer, intent (in) :: iv, imu, is

    select case (vms_layout)
    case ('vms')
       idx_vmu = iv+lo%nvgrid + lo%nvpa*(imu-1 + lo%nmu*(is-1))
    case ('mvs')
       idx_vmu = imu-1 + lo%nmu*(iv+lo%nvgrid + lo%nvpa*(is-1))
    end select

  end function idx_vmu

  elemental function idx_local_vmu (lo, iv, imu, is)

    implicit none
    logical :: idx_local_vmu
    type (vmu_layout_type), intent (in) :: lo
    integer, intent (in) :: iv, imu, is

    idx_local_vmu = idx_local(lo, idx(lo, iv, imu, is))
  end function idx_local_vmu

  elemental function ig_local_vmu (lo, ig)
    implicit none
    logical :: ig_local_vmu
    type (vmu_layout_type), intent (in) :: lo
    integer, intent (in) :: ig

    ig_local_vmu = lo%iproc == proc_id(lo, ig)
  end function ig_local_vmu

  elemental subroutine kxkyzidx2vmuidx (iv, imu, ikxkyz, kxkyz_lo, vmu_lo, iky, ikx, ig, ivmu)
    implicit none
    integer, intent (in) :: iv, imu, ikxkyz
    type (kxkyz_layout_type), intent (in) :: kxkyz_lo
    type (vmu_layout_type), intent (in) :: vmu_lo
    integer, intent (out) :: iky, ikx, ig, ivmu

    iky = iky_idx(kxkyz_lo,ikxkyz)
    ikx = ikx_idx(kxkyz_lo,ikxkyz)
    ig = iz_idx(kxkyz_lo,ikxkyz)
    ivmu = idx(vmu_lo, iv, imu, is_idx(kxkyz_lo,ikxkyz))
  end subroutine kxkyzidx2vmuidx

  elemental subroutine kxyzidx2vmuidx (iv, imu, ikxyz, kxyz_lo, vmu_lo, iy, ikx, ig, ivmu)
    implicit none
    integer, intent (in) :: iv, imu, ikxyz
    type (kxyz_layout_type), intent (in) :: kxyz_lo
    type (vmu_layout_type), intent (in) :: vmu_lo
    integer, intent (out) :: iy, ikx, ig, ivmu

    iy = iy_idx(kxyz_lo,ikxyz)
    ikx = ikx_idx(kxyz_lo,ikxyz)
    ig = iz_idx(kxyz_lo,ikxyz)
    ivmu = idx(vmu_lo, iv, imu, is_idx(kxyz_lo,ikxyz))
  end subroutine kxyzidx2vmuidx

  elemental subroutine xyzidx2vmuidx (iv, imu, ixyz, xyz_lo, vmu_lo, iy, ix, iz, ivmu)
    implicit none
    integer, intent (in) :: iv, imu, ixyz
    type (xyz_layout_type), intent (in) :: xyz_lo
    type (vmu_layout_type), intent (in) :: vmu_lo
    integer, intent (out) :: iy, ix, iz, ivmu

    iy = iy_idx(xyz_lo,ixyz)
    ix = ix_idx(xyz_lo,ixyz)
    iz = iz_idx(xyz_lo,ixyz)
    ivmu = idx(vmu_lo, iv, imu, is_idx(xyz_lo,ixyz))
  end subroutine xyzidx2vmuidx

  subroutine finish_layouts
    
    implicit none
    
  end subroutine finish_layouts
  
end module stella_layouts
