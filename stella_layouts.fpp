# include "define.inc"

module stella_layouts

  use layouts_type, only: gvmu_layout_type, gxyz_layout_type, gy_layout_type
  
  implicit none

  private
  
  public :: layout, finish_layouts
  
  public :: init_stella_layouts, init_dist_fn_layouts
  public :: gvmu_lo, gy_lo, gxyz_lo
  
  public :: init_x_transform_layouts, init_y_transform_layouts
  public :: xxf_lo, xxf_layout_type, yxf_lo, yxf_layout_type
  public :: vmuidx2yidx, vmuidx2zkxkyidx
  public :: xxfidx2yxfidx, yxfidx2xxfidx
  public :: xxf_ky_is_zero
  
  public :: ig_idx, iky_idx, ikx_idx, iv_idx, imu_idx, is_idx
  public :: idx, proc_id, idx_local
  
  logical :: initialized_x_transform = .false.
  logical :: initialized_y_transform = .false.
  
  character (len=6) :: layout
  logical :: exist
  
  type (gvmu_layout_type) :: gvmu_lo
  type (gy_layout_type) :: gy_lo
  type (gxyz_layout_type) :: gxyz_lo
  
  type :: xxf_layout_type
     integer :: iproc
     integer :: ntgrid, nvgrid, naky, ntheta0, nx, nadd, nmu, nspec, ntgridtotal
     integer :: llim_world, ulim_world, llim_proc, ulim_proc, ulim_alloc, blocksize, gsize
     integer :: llim_group, ulim_group, igroup, ngroup, nprocset, iset, nset, groupblocksize
     integer :: small_block_size, block_multiple, large_block_size, num_small, num_large
     integer :: small_block_balance_factor, large_block_balance_factor
  end type xxf_layout_type

  type (xxf_layout_type) :: xxf_lo

  type :: yxf_layout_type
     integer :: iproc
     integer :: ntgrid, nvgrid, naky, ny, ntheta0, nx, nmu, nspec, ntgridtotal, nvgridtotal
     integer :: llim_world, ulim_world, llim_proc, ulim_proc, ulim_alloc, blocksize, gsize
     integer :: llim_group, ulim_group, igroup, ngroup, nprocset, iset, nset, groupblocksize
     integer :: small_block_size, block_multiple, large_block_size, num_small, num_large
     integer :: small_block_balance_factor, large_block_balance_factor
     integer :: ikx_ord, ig_ord, iv_ord, imu_ord, is_ord
     integer :: ikx_comp, ig_comp, iv_comp, imu_comp, is_comp
     integer, dimension (5) :: compound_count
  end type yxf_layout_type

  type (yxf_layout_type) :: yxf_lo

  interface ig_idx
     module procedure ig_idx_gvmu
     module procedure ig_idx_gy
     module procedure ig_idx_xxf
     module procedure ig_idx_yxf
  end interface

  interface iv_idx
     module procedure iv_idx_gy
     module procedure iv_idx_gxyz
     module procedure iv_idx_xxf
     module procedure iv_idx_yxf
  end interface

  interface iky_idx
     module procedure iky_idx_gvmu
     module procedure iky_idx_xxf
  end interface

  interface ikx_idx
     module procedure ikx_idx_gy
     module procedure ikx_idx_gvmu
     module procedure ikx_idx_yxf
  end interface

  interface imu_idx
     module procedure imu_idx_gy
     module procedure imu_idx_gxyz
     module procedure imu_idx_xxf
     module procedure imu_idx_yxf
  end interface

  interface is_idx
     module procedure is_idx_gvmu
     module procedure is_idx_gy
     module procedure is_idx_gxyz
     module procedure is_idx_xxf
     module procedure is_idx_yxf
  end interface

  interface proc_id
     module procedure proc_id_gy
     module procedure proc_id_gvmu
     module procedure proc_id_gxyz
     module procedure proc_id_xxf
     module procedure proc_id_yxf
  end interface

  interface idx
     module procedure idx_gy
     module procedure idx_gvmu
     module procedure idx_gxyz
     module procedure idx_xxf
     module procedure idx_yxf
  end interface

  interface idx_local
     module procedure idx_local_gy, ig_local_gy
     module procedure idx_local_gvmu, ig_local_gvmu
     module procedure idx_local_gxyz, ig_local_gxyz
     module procedure idx_local_xxf, ig_local_xxf
     module procedure idx_local_yxf, ig_local_yxf
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
    use file_utils, only: input_unit, error_unit, input_unit_exist, error_unit
    implicit none
    integer :: in_file
    namelist /layouts_knobs/ layout
    layout = 'vmyxzs'
    in_file=input_unit_exist("layouts_knobs", exist)
    if (exist) read (unit=input_unit("layouts_knobs"), nml=layouts_knobs)
    if (layout.ne.'xyzvms' .and. &
         layout.ne.'vmyxzs' .and. &
         layout.ne.'xyzmvs' .and. &
         layout.ne.'zyxvms') then
       write(6,*) "stella_layouts: read_parameters finds illegal layout=",layout," =>stop"
       stop
    endif

  end subroutine read_parameters
    
  subroutine broadcast_results
    use mp, only: broadcast
    implicit none

    call broadcast (layout)

  end subroutine broadcast_results

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Distribution function layouts
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine init_dist_fn_layouts (ntgrid, naky, nakx, nvgrid, nmu, nspec, ny)

    implicit none

    integer, intent (in) :: ntgrid, naky, nakx, nvgrid, nmu, nspec, ny

    call init_gvmu_layout (ntgrid, naky, nakx, nvgrid, nmu, nspec)
    call init_gxyz_layout (ntgrid, naky, nakx, nvgrid, nmu, nspec)
    call init_gy_layout (ntgrid, naky, nakx, nvgrid, nmu, nspec, ny)

  end subroutine init_dist_fn_layouts

  subroutine init_gvmu_layout &
       (ntgrid, naky, nakx, nvgrid, nmu, nspec)
    
    use mp, only: iproc, nproc

    implicit none

    integer, intent (in) :: ntgrid, naky, nakx, nvgrid, nmu, nspec
    logical, save :: initialized = .false.

    if (initialized) return
    initialized = .true.
   
    gvmu_lo%iproc = iproc
    gvmu_lo%ntgrid = ntgrid
    gvmu_lo%ntheta = 2*ntgrid+1
    gvmu_lo%naky = naky
    gvmu_lo%nakx = nakx
    gvmu_lo%nvgrid = nvgrid
    gvmu_lo%nvpa = 2*nvgrid+1
    gvmu_lo%nmu = nmu
    gvmu_lo%nspec = nspec
    gvmu_lo%llim_world = 0
    gvmu_lo%ulim_world = naky*nakx*gvmu_lo%ntheta*nspec - 1
    gvmu_lo%blocksize = gvmu_lo%ulim_world/nproc + 1
    gvmu_lo%llim_proc = gvmu_lo%blocksize*iproc
    gvmu_lo%ulim_proc = min(gvmu_lo%ulim_world, gvmu_lo%llim_proc + gvmu_lo%blocksize - 1)
    gvmu_lo%ulim_alloc = max(gvmu_lo%llim_proc, gvmu_lo%ulim_proc)

  end subroutine init_gvmu_layout

  elemental function is_idx_gvmu (lo, i)
    implicit none
    integer :: is_idx_gvmu
    type (gvmu_layout_type), intent (in) :: lo
    integer, intent (in) :: i
    is_idx_gvmu = 1 + mod((i - lo%llim_world)/lo%nakx/lo%naky/lo%ntheta, lo%nspec)
  end function is_idx_gvmu

  elemental function ikx_idx_gvmu (lo, i)

    implicit none

    integer :: ikx_idx_gvmu
    type (gvmu_layout_type), intent (in) :: lo
    integer, intent (in) :: i

    select case (layout)
    case ('vmyxzs')
       ikx_idx_gvmu = 1 + mod((i - lo%llim_world)/lo%naky, lo%nakx)
    case ('xyzvms')
       ikx_idx_gvmu = 1 + mod((i - lo%llim_world), lo%nakx)
    case ('xyzmvs')
       ikx_idx_gvmu = 1 + mod((i - lo%llim_world), lo%nakx)
    case ('zyxvms')
       ikx_idx_gvmu = 1 + mod((i - lo%llim_world)/lo%naky/lo%ntheta, lo%nakx)
    end select

  end function ikx_idx_gvmu

  elemental function iky_idx_gvmu (lo, i)

    implicit none
    integer :: iky_idx_gvmu
    type (gvmu_layout_type), intent (in) :: lo
    integer, intent (in) :: i

    select case (layout)
    case ('vmyxzs')
       iky_idx_gvmu = 1 + mod(i - lo%llim_world, lo%naky)
    case ('xyzvms')
       iky_idx_gvmu = 1 + mod((i - lo%llim_world)/lo%nakx, lo%naky)
    case ('xyzmvs')
       iky_idx_gvmu = 1 + mod((i - lo%llim_world)/lo%nakx, lo%naky)
    case ('zyxvms')
       iky_idx_gvmu = 1 + mod((i - lo%llim_world)/lo%ntheta, lo%naky)
    end select

  end function iky_idx_gvmu

  elemental function ig_idx_gvmu (lo, i)

    implicit none
    integer :: ig_idx_gvmu
    type (gvmu_layout_type), intent (in) :: lo
    integer, intent (in) :: i

    select case (layout)
    case ('vmyxzs')
       ig_idx_gvmu = -lo%ntgrid + mod((i - lo%llim_world)/lo%nakx/lo%naky, lo%ntheta)
    case ('xyzvms')
       ig_idx_gvmu = -lo%ntgrid + mod((i - lo%llim_world)/lo%nakx/lo%naky, lo%ntheta)
    case ('xyzmvs')
       ig_idx_gvmu = -lo%ntgrid + mod((i - lo%llim_world)/lo%nakx/lo%naky, lo%ntheta)
    case ('zyxvms')
       ig_idx_gvmu = -lo%ntgrid + mod((i - lo%llim_world), lo%ntheta)
    end select

  end function ig_idx_gvmu

  elemental function proc_id_gvmu (lo, i)
    implicit none
    integer :: proc_id_gvmu
    type (gvmu_layout_type), intent (in) :: lo
    integer, intent (in) :: i

    proc_id_gvmu = i/lo%blocksize

  end function proc_id_gvmu

  elemental function idx_gvmu (lo, iky, ikx, ig, is)

    implicit none

    integer :: idx_gvmu
    type (gvmu_layout_type), intent (in) :: lo
    integer, intent (in) :: iky, ikx, ig, is

    select case (layout)
    case ('vmyxzs')
       idx_gvmu = iky-1 + lo%naky*(ikx-1 + lo%nakx*(ig+lo%ntgrid + lo%ntheta*(is-1)))
    case ('xyzvms')
       idx_gvmu = ikx-1 + lo%nakx*(iky-1 + lo%naky*(ig+lo%ntgrid + lo%ntheta*(is-1)))
    case ('xyzmvs')
       idx_gvmu = ikx-1 + lo%nakx*(iky-1 + lo%naky*(ig+lo%ntgrid + lo%ntheta*(is-1)))
    case ('zyxvms')
       idx_gvmu = ig+lo%ntgrid + lo%ntheta*(iky-1 + lo%naky*(ikx-1))
    end select

  end function idx_gvmu

  elemental function idx_local_gvmu (lo, iky, ikx, ig, is)

    implicit none
    logical :: idx_local_gvmu
    type (gvmu_layout_type), intent (in) :: lo
    integer, intent (in) :: iky, ikx, ig, is

    idx_local_gvmu = idx_local(lo, idx(lo, iky, ikx, ig, is))
  end function idx_local_gvmu

  elemental function ig_local_gvmu (lo, ig)
    implicit none
    logical :: ig_local_gvmu
    type (gvmu_layout_type), intent (in) :: lo
    integer, intent (in) :: ig

    ig_local_gvmu = lo%iproc == proc_id(lo, ig)
  end function ig_local_gvmu

  subroutine init_gy_layout &
       (ntgrid, naky, nakx, nvgrid, nmu, nspec, ny)
    
    use mp, only: iproc, nproc

    implicit none

    integer, intent (in) :: ntgrid, naky, nakx, nvgrid, nmu, nspec, ny
    logical, save :: initialized = .false.

    if (initialized) return
    initialized = .true.
   
    gy_lo%iproc = iproc
    gy_lo%ntheta = 2*ntgrid+1
    gy_lo%ntgrid = ntgrid
    gy_lo%naky = naky
    gy_lo%nakx = nakx
    gy_lo%nvgrid = nvgrid
    gy_lo%nvpa = 2*nvgrid+1
    gy_lo%nmu = nmu
    gy_lo%nspec = nspec
    gy_lo%ny = ny
    gy_lo%llim_world = 0
    gy_lo%ulim_world = gy_lo%nvpa*gy_lo%ntheta*nakx*nmu*nspec - 1

    gy_lo%ikx_ord=1
    gy_lo%ig_ord=2
    gy_lo%iv_ord=3
    gy_lo%imu_ord=4
    gy_lo%is_ord=5

    gy_lo%compound_count(1)=1
    gy_lo%compound_count(2)=gy_lo%nakx
    gy_lo%compound_count(3)=gy_lo%compound_count(2)*gy_lo%ntheta
    gy_lo%compound_count(4)=gy_lo%compound_count(3)*gy_lo%nvpa
    gy_lo%compound_count(5)=gy_lo%compound_count(4)*nmu

    gy_lo%ig_comp=gy_lo%compound_count(gy_lo%ig_ord)
    gy_lo%iv_comp=gy_lo%compound_count(gy_lo%iv_ord)
    gy_lo%ikx_comp=gy_lo%compound_count(gy_lo%ikx_ord)
    gy_lo%imu_comp=gy_lo%compound_count(gy_lo%imu_ord)
    gy_lo%is_comp=gy_lo%compound_count(gy_lo%is_ord)
      
    gy_lo%blocksize = gy_lo%ulim_world/nproc + 1
    gy_lo%llim_proc = gy_lo%blocksize*iproc
    gy_lo%ulim_proc = min(gy_lo%ulim_world, gy_lo%llim_proc + gy_lo%blocksize - 1)
    gy_lo%ulim_alloc = max(gy_lo%llim_proc, gy_lo%ulim_proc)

  end subroutine init_gy_layout

  elemental function is_idx_gy (lo, i)
    implicit none
    integer :: is_idx_gy
    type (gy_layout_type), intent (in) :: lo
    integer, intent (in) :: i
    is_idx_gy = 1 + mod((i-lo%llim_world)/lo%is_comp, lo%nspec)
  end function is_idx_gy

  elemental function imu_idx_gy (lo, i)
    implicit none
    integer :: imu_idx_gy
    type (gy_layout_type), intent (in) :: lo
    integer, intent (in) :: i
    imu_idx_gy = 1 + mod((i-lo%llim_world)/lo%imu_comp, lo%nmu)
  end function imu_idx_gy

  elemental function iv_idx_gy (lo, i)
    implicit none
    integer :: iv_idx_gy
    type (gy_layout_type), intent (in) :: lo
    integer, intent (in) :: i
    iv_idx_gy = -lo%nvgrid + mod((i-lo%llim_world)/lo%iv_comp, lo%nvpa)
  end function iv_idx_gy

  elemental function ig_idx_gy (lo, i)
    implicit none
    integer :: ig_idx_gy
    type (gy_layout_type), intent (in) :: lo
    integer, intent (in) :: i
    ig_idx_gy = -lo%ntgrid + mod((i-lo%llim_world)/lo%ig_comp, lo%ntheta)
  end function ig_idx_gy

  elemental function ikx_idx_gy (lo, i)
    implicit none
    integer :: ikx_idx_gy
    type (gy_layout_type), intent (in) :: lo
    integer, intent (in) :: i
    ikx_idx_gy = 1 + mod((i-lo%llim_world)/lo%ikx_comp, lo%nakx)
  end function ikx_idx_gy

  elemental function idx_gy (lo, ig, iv, ikx, imu, is)
    implicit none
    integer :: idx_gy
    type (gy_layout_type), intent (in) :: lo
    integer, intent (in) :: ig, iv, ikx, imu, is
    idx_gy = ikx-1 + lo%nakx*(ig+lo%ntgrid + lo%ntheta*(iv+lo%nvgrid &
         + lo%nvpa*(imu-1 + lo%nmu*(is-1))))
  end function idx_gy

  elemental function proc_id_gy (lo, i)
    implicit none
    integer :: proc_id_gy
    type (gy_layout_type), intent (in) :: lo
    integer, intent (in) :: i
    proc_id_gy = i/lo%blocksize
  end function proc_id_gy

  elemental function idx_local_gy (lo, ig, iv, ikx, imu, is)
    implicit none
    logical :: idx_local_gy
    type (gy_layout_type), intent (in) :: lo
    integer, intent (in) :: ig, iv, ikx, imu, is
    idx_local_gy = idx_local (lo, idx(lo, ig, iv, ikx, imu, is))
  end function idx_local_gy

  elemental function ig_local_gy (lo, i)
    implicit none
    logical ig_local_gy
    type (gy_layout_type), intent (in) :: lo
    integer, intent (in) :: i
    ig_local_gy = lo%iproc == proc_id(lo, i)
  end function ig_local_gy

  subroutine init_gxyz_layout &
       (ntgrid, naky, nakx, nvgrid, nmu, nspec)
    
    use mp, only: iproc, nproc

    implicit none

    integer, intent (in) :: ntgrid, naky, nakx, nvgrid, nmu, nspec
    logical, save :: initialized = .false.

    if (initialized) return
    initialized = .true.
   
    gxyz_lo%xyz = .true.
    gxyz_lo%iproc = iproc
    gxyz_lo%ntheta = 2*ntgrid+1
    gxyz_lo%ntgrid = ntgrid
    gxyz_lo%naky = naky
    gxyz_lo%nakx = nakx
    gxyz_lo%nvgrid = nvgrid
    gxyz_lo%nvpa = 2*nvgrid+1
    gxyz_lo%nmu = nmu
    gxyz_lo%nspec = nspec
    gxyz_lo%llim_world = 0
    gxyz_lo%ulim_world = gxyz_lo%nvpa*nmu*nspec - 1
      
    gxyz_lo%blocksize = gxyz_lo%ulim_world/nproc + 1
    gxyz_lo%llim_proc = gxyz_lo%blocksize*iproc
    gxyz_lo%ulim_proc = min(gxyz_lo%ulim_world, gxyz_lo%llim_proc + gxyz_lo%blocksize - 1)
    gxyz_lo%ulim_alloc = max(gxyz_lo%llim_proc, gxyz_lo%ulim_proc)

  end subroutine init_gxyz_layout

  elemental function is_idx_gxyz (lo, i)

    implicit none
    integer :: is_idx_gxyz
    type (gxyz_layout_type), intent (in) :: lo
    integer, intent (in) :: i

    ! the order of the division does not matter, so no need for branching
    is_idx_gxyz = 1 + mod((i - lo%llim_world)/lo%nvpa/lo%nmu, lo%nspec)

  end function is_idx_gxyz

  elemental function imu_idx_gxyz (lo, i)

    implicit none

    integer :: imu_idx_gxyz
    type (gxyz_layout_type), intent (in) :: lo
    integer, intent (in) :: i

    select case (layout)
    case ('vmyxzs')
       imu_idx_gxyz = 1 + mod((i - lo%llim_world)/lo%nvpa, lo%nmu)
    case ('xyzvms')
       imu_idx_gxyz = 1 + mod((i - lo%llim_world)/lo%nvpa, lo%nmu)
    case ('xyzmvs')
       imu_idx_gxyz = 1 + mod((i - lo%llim_world), lo%nmu)
    case ('zyxvms')
       imu_idx_gxyz = 1 + mod((i - lo%llim_world)/lo%nvpa, lo%nmu)
    end select

  end function imu_idx_gxyz

  elemental function iv_idx_gxyz (lo, i)

    implicit none
    integer :: iv_idx_gxyz
    type (gxyz_layout_type), intent (in) :: lo
    integer, intent (in) :: i

    select case (layout)
    case ('vmyxzs')
       iv_idx_gxyz = -lo%nvgrid + mod((i - lo%llim_world), lo%nmu)
    case ('xyzvms')
       iv_idx_gxyz = -lo%nvgrid + mod((i - lo%llim_world), lo%nvpa)
    case ('xyzmvs')
       iv_idx_gxyz = -lo%nvgrid + mod((i - lo%llim_world)/lo%nmu, lo%nvpa)
    case ('zyxvms')
       iv_idx_gxyz = -lo%nvgrid + mod((i - lo%llim_world), lo%nvpa)
    end select

  end function iv_idx_gxyz

  elemental function proc_id_gxyz (lo, i)
    implicit none
    integer :: proc_id_gxyz
    type (gxyz_layout_type), intent (in) :: lo
    integer, intent (in) :: i

    proc_id_gxyz = i/lo%blocksize

  end function proc_id_gxyz

  elemental function idx_gxyz (lo, iv, imu, is)

    implicit none

    integer :: idx_gxyz
    type (gxyz_layout_type), intent (in) :: lo
    integer, intent (in) :: iv, imu, is

    select case (layout)
    case ('vmyxzs')
       idx_gxyz = iv+lo%nvgrid + lo%nvpa*(imu-1 + lo%nmu*(is-1))
    case ('xyzvms')
       idx_gxyz = iv+lo%nvgrid + lo%nvpa*(imu-1 + lo%nmu*(is-1))
    case ('xyzmvs')
       idx_gxyz = imu-1 + lo%nmu*(iv+lo%nvgrid + lo%nvpa*(is-1))
    case ('zyxvms')
       idx_gxyz = iv+lo%nvgrid + lo%nvpa*(imu-1 + lo%nmu*(is-1))
    end select

  end function idx_gxyz

  elemental function idx_local_gxyz (lo, iv, imu, is)

    implicit none
    logical :: idx_local_gxyz
    type (gxyz_layout_type), intent (in) :: lo
    integer, intent (in) :: iv, imu, is

    idx_local_gxyz = idx_local(lo, idx(lo, iv, imu, is))
  end function idx_local_gxyz

  elemental function ig_local_gxyz (lo, ig)
    implicit none
    logical :: ig_local_gxyz
    type (gxyz_layout_type), intent (in) :: lo
    integer, intent (in) :: ig

    ig_local_gxyz = lo%iproc == proc_id(lo, ig)
  end function ig_local_gxyz

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! X-space layouts
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  subroutine init_x_transform_layouts &
       (ntgrid, naky, ntheta0, nvgrid, nmu, nspec, nx)
    use mp, only: iproc, nproc, barrier
    implicit none
    integer, intent (in) :: ntgrid, nvgrid, naky, ntheta0, nmu, nspec, nx

    if (initialized_x_transform) return
    initialized_x_transform = .true.

    xxf_lo%iproc = iproc
    xxf_lo%ntgrid = ntgrid
    xxf_lo%ntgridtotal = (2*ntgrid+1)
    xxf_lo%nvgrid = nvgrid
    xxf_lo%naky = naky
    xxf_lo%ntheta0 = ntheta0
    if (nx > ntheta0) then
       xxf_lo%nx = nx
    else
       xxf_lo%nx = (3*ntheta0+1)/2
    end if
    xxf_lo%nadd = xxf_lo%nx - ntheta0
    xxf_lo%nmu = nmu
    xxf_lo%nspec = nspec
    xxf_lo%llim_world = 0
    xxf_lo%ulim_world = naky*(2*ntgrid+1)*(2*nvgrid+1)*nmu*nspec - 1

    xxf_lo%blocksize = xxf_lo%ulim_world/nproc + 1
    xxf_lo%llim_proc = xxf_lo%blocksize*iproc
    xxf_lo%ulim_proc &
         = min(xxf_lo%ulim_world, xxf_lo%llim_proc + xxf_lo%blocksize - 1)
    xxf_lo%ulim_alloc = max(xxf_lo%llim_proc, xxf_lo%ulim_proc)

  end subroutine init_x_transform_layouts

  elemental function is_idx_xxf (lo, i)
    implicit none
    integer :: is_idx_xxf
    type (xxf_layout_type), intent (in) :: lo
    integer, intent (in) :: i

    ! TT: the order of the division does not matter, so no need for branching
    is_idx_xxf = 1 + mod((i - lo%llim_world)/lo%naky/(2*lo%ntgrid + 1)/(2*lo%nvgrid+1)/lo%nmu, lo%nspec)

  end function is_idx_xxf

  elemental function imu_idx_xxf (lo, i)
    implicit none
    integer :: imu_idx_xxf
    type (xxf_layout_type), intent (in) :: lo
    integer, intent (in) :: i

    imu_idx_xxf = 1 + mod((i - lo%llim_world)/lo%naky/(2*lo%ntgrid + 1)/(2*lo%nvgrid + 1), lo%nmu)
  end function imu_idx_xxf

  elemental function iv_idx_xxf (lo, i)
    implicit none
    integer :: iv_idx_xxf
    type (xxf_layout_type), intent (in) :: lo
    integer, intent (in) :: i

    iv_idx_xxf = -lo%nvgrid + mod((i - lo%llim_world)/lo%naky/(2*lo%ntgrid + 1), 2*lo%nvgrid + 1)
  end function iv_idx_xxf

  elemental function ig_idx_xxf (lo, i)
    implicit none
    integer :: ig_idx_xxf
    type (xxf_layout_type), intent (in) :: lo
    integer, intent (in) :: i

    ig_idx_xxf = -lo%ntgrid + mod((i - lo%llim_world)/lo%naky, 2*lo%ntgrid + 1)
  end function ig_idx_xxf

  elemental function iky_idx_xxf (lo, i)
    implicit none
    integer :: iky_idx_xxf
    type (xxf_layout_type), intent (in) :: lo
    integer, intent (in) :: i

    iky_idx_xxf = 1 + mod(i - lo%llim_world, lo%naky)
  end function iky_idx_xxf

  elemental function idx_xxf (lo, ig, iv, iky, imu, is)
    implicit none
    integer :: idx_xxf
    type (xxf_layout_type), intent (in) :: lo
    integer, intent (in) :: ig, iv, iky, imu, is

    idx_xxf = iky-1 + lo%naky*(ig+lo%ntgrid + (2*lo%ntgrid+1)*(iv+lo%nvgrid &
         + (2*lo%nvgrid+1)*(imu-1 + lo%nmu*(is-1))))
  end function idx_xxf

  elemental function proc_id_xxf (lo, i)
    implicit none
    integer :: proc_id_xxf
    type (xxf_layout_type), intent (in) :: lo
    integer, intent (in) :: i
    
    proc_id_xxf = i/lo%blocksize

  end function proc_id_xxf

  elemental function idx_local_xxf (lo, ig, iv, iky, imu, is)
    implicit none
    logical :: idx_local_xxf
    type (xxf_layout_type), intent (in) :: lo
    integer, intent (in) :: ig, iv, iky, imu, is
    idx_local_xxf = idx_local (lo, idx(lo, ig, iv, iky, imu, is))
  end function idx_local_xxf

  elemental function ig_local_xxf (lo, i)
    implicit none
    logical ig_local_xxf
    type (xxf_layout_type), intent (in) :: lo
    integer, intent (in) :: i
    ig_local_xxf = lo%iproc == proc_id(lo, i)
  end function ig_local_xxf

  elemental function xxf_ky_is_zero (lo, i)
    implicit none
    logical xxf_ky_is_zero
    type (xxf_layout_type), intent (in) :: lo
    integer, intent (in) :: i

    xxf_ky_is_zero = 0 == mod(i, lo%naky)

  end function xxf_ky_is_zero

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Y-space layouts
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine init_y_transform_layouts &
       (ntgrid, naky, ntheta0, nvgrid, nmu, nspec, nx, ny)
    use mp, only: iproc, nproc
    implicit none
    integer, intent (in) :: ntgrid, nvgrid, naky, ntheta0, nmu, nspec
    integer, intent (in) :: nx, ny
    integer :: nnx, nny

    if (initialized_y_transform) return
    initialized_y_transform = .true.

    if (nx > ntheta0) then
       nnx = nx
    else
       nnx = (3*ntheta0+1)/2
    end if
    if (ny > naky) then
       nny = ny
    else
       nny = 3*naky
    end if

    yxf_lo%iproc = iproc
    yxf_lo%ntgrid = ntgrid
    yxf_lo%ntgridtotal = (2*ntgrid+1)
    yxf_lo%nvgrid = nvgrid
    yxf_lo%nvgridtotal = (2*nvgrid+1)
    yxf_lo%naky = naky
    yxf_lo%ny = nny
    yxf_lo%ntheta0 = ntheta0
    yxf_lo%nx = nnx
    yxf_lo%nmu = nmu
    yxf_lo%nspec = nspec
    yxf_lo%llim_world = 0
    yxf_lo%ulim_world = nnx*(2*ntgrid+1)*(2*nvgrid+1)*nmu*nspec - 1

    yxf_lo%ikx_ord=1
    yxf_lo%ig_ord=2
    yxf_lo%iv_ord=3
    yxf_lo%imu_ord=4
    yxf_lo%is_ord=5

    yxf_lo%compound_count(1)=1
    yxf_lo%compound_count(2)=yxf_lo%nx
    yxf_lo%compound_count(3)=yxf_lo%compound_count(2)*yxf_lo%ntgridtotal
    yxf_lo%compound_count(4)=yxf_lo%compound_count(3)*yxf_lo%nvgridtotal
    yxf_lo%compound_count(5)=yxf_lo%compound_count(4)*nmu

    yxf_lo%ig_comp=yxf_lo%compound_count(yxf_lo%ig_ord)
    yxf_lo%iv_comp=yxf_lo%compound_count(yxf_lo%iv_ord)
    yxf_lo%ikx_comp=yxf_lo%compound_count(yxf_lo%ikx_ord)
    yxf_lo%imu_comp=yxf_lo%compound_count(yxf_lo%imu_ord)
    yxf_lo%is_comp=yxf_lo%compound_count(yxf_lo%is_ord)

    yxf_lo%blocksize = yxf_lo%ulim_world/nproc + 1
    yxf_lo%llim_proc = yxf_lo%blocksize*iproc
    yxf_lo%ulim_proc &
         = min(yxf_lo%ulim_world, yxf_lo%llim_proc + yxf_lo%blocksize - 1)
    yxf_lo%ulim_alloc = max(yxf_lo%llim_proc, yxf_lo%ulim_proc)
    
  end subroutine init_y_transform_layouts

  elemental function is_idx_yxf (lo, i)
    implicit none
    integer :: is_idx_yxf
    type (yxf_layout_type), intent (in) :: lo
    integer, intent (in) :: i
    is_idx_yxf = 1 + mod((i-lo%llim_world)/lo%is_comp, lo%nspec)
  end function is_idx_yxf

  elemental function imu_idx_yxf (lo, i)
    implicit none
    integer :: imu_idx_yxf
    type (yxf_layout_type), intent (in) :: lo
    integer, intent (in) :: i
    imu_idx_yxf = 1 + mod((i-lo%llim_world)/lo%imu_comp, lo%nmu)
  end function imu_idx_yxf

  elemental function iv_idx_yxf (lo, i)
    implicit none
    integer :: iv_idx_yxf
    type (yxf_layout_type), intent (in) :: lo
    integer, intent (in) :: i
    iv_idx_yxf = -lo%nvgrid + mod((i-lo%llim_world)/lo%iv_comp, lo%nvgridtotal)
  end function iv_idx_yxf

  elemental function ig_idx_yxf (lo, i)
    implicit none
    integer :: ig_idx_yxf
    type (yxf_layout_type), intent (in) :: lo
    integer, intent (in) :: i
    ig_idx_yxf = -lo%ntgrid + mod((i-lo%llim_world)/lo%ig_comp, lo%ntgridtotal)
  end function ig_idx_yxf

  elemental function ikx_idx_yxf (lo, i)
    implicit none
    integer :: ikx_idx_yxf
    type (yxf_layout_type), intent (in) :: lo
    integer, intent (in) :: i
    ikx_idx_yxf = 1 + mod((i-lo%llim_world)/lo%ikx_comp, lo%nx)
  end function ikx_idx_yxf

  elemental function idx_yxf (lo, ig, iv, ikx, imu, is)
    implicit none
    integer :: idx_yxf
    type (yxf_layout_type), intent (in) :: lo
    integer, intent (in) :: ig, iv, ikx, imu, is

    idx_yxf = ikx-1 + lo%nx*(ig+lo%ntgrid + (2*lo%ntgrid+1)*(iv+lo%nvgrid &
         + (2*lo%nvgrid+1)*(imu-1 + lo%nmu*(is-1))))
  end function idx_yxf

  elemental function proc_id_yxf (lo, i)
    implicit none
    integer :: proc_id_yxf
    type (yxf_layout_type), intent (in) :: lo
    integer, intent (in) :: i

    proc_id_yxf = i/lo%blocksize

  end function proc_id_yxf

  elemental function idx_local_yxf (lo, ig, iv, ikx, imu, is)
    implicit none
    logical :: idx_local_yxf
    type (yxf_layout_type), intent (in) :: lo
    integer, intent (in) :: ig, iv, ikx, imu, is
    idx_local_yxf = idx_local (lo, idx(lo, ig, iv, ikx, imu, is))
  end function idx_local_yxf

  elemental function ig_local_yxf (lo, i)
    implicit none
    logical ig_local_yxf
    type (yxf_layout_type), intent (in) :: lo
    integer, intent (in) :: i
    ig_local_yxf = lo%iproc == proc_id(lo, i)
  end function ig_local_yxf

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Transformation subroutines
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


  ! DO NOT USE WITHOUT FIXING TO ACCOUNT FOR G_LO TO GXYZ_LO
!   elemental subroutine gidx2xxfidx (ig, iv, it_in, iglo, gxyz_lo, xxf_lo, it_out, ixxf)

!     implicit none
!     integer, intent (in) :: ig, iv, it_in, iglo
!     type (gxyz_layout_type), intent (in) :: gxyz_lo
!     type (xxf_layout_type), intent (in) :: xxf_lo
!     integer, intent (out) :: ixxf, it_out

!     if (it_in > (xxf_lo%ntheta0+1)/2) then
!        it_out = it_in - xxf_lo%ntheta0 + xxf_lo%nx
!     else
!        it_out = it_in
!     end if

!     ixxf = idx(xxf_lo, ig, iv, iky_idx(gxyz_lo,iglo), &
!          imu_idx(gxyz_lo,iglo), is_idx(gxyz_lo,iglo))
!   end subroutine gidx2xxfidx

  elemental subroutine vmuidx2yidx (iv, imu, ivmu, gvmu_lo, gy_lo, iky, iy)
    implicit none
    integer, intent (in) :: iv, imu, ivmu
    type (gvmu_layout_type), intent (in) :: gvmu_lo
    type (gy_layout_type), intent (in) :: gy_lo
    integer, intent (out) :: iky, iy

    iky = iky_idx(gvmu_lo,ivmu)
    iy = idx(gy_lo, ig_idx(gvmu_lo,ivmu), iv, &
         ikx_idx(gvmu_lo,ivmu), imu, is_idx(gvmu_lo,ivmu))
  end subroutine vmuidx2yidx

  elemental subroutine vmuidx2zkxkyidx (iv, imu, ivmu, gvmu_lo, gxyz_lo, iky, ikx, ig, ixyz)
    implicit none
    integer, intent (in) :: iv, imu, ivmu
    type (gvmu_layout_type), intent (in) :: gvmu_lo
    type (gxyz_layout_type), intent (in) :: gxyz_lo
    integer, intent (out) :: iky, ikx, ig, ixyz

    iky = iky_idx(gvmu_lo,ivmu)
    ikx = ikx_idx(gvmu_lo,ivmu)
    ig = ig_idx(gvmu_lo,ivmu)
    ixyz = idx(gxyz_lo, iv, imu, is_idx(gvmu_lo,ivmu))
  end subroutine vmuidx2zkxkyidx

  elemental subroutine xxfidx2yxfidx (ikx, ixxf, xxf_lo, yxf_lo, iky, iyxf)
    implicit none
    integer, intent (in) :: ikx, ixxf
    type (xxf_layout_type), intent (in) :: xxf_lo
    type (yxf_layout_type), intent (in) :: yxf_lo
    integer, intent (out) :: iky, iyxf

    iky = iky_idx(xxf_lo,ixxf)
    iyxf = idx(yxf_lo, ig_idx(xxf_lo,ixxf), iv_idx(xxf_lo,ixxf), &
         ikx, imu_idx(xxf_lo,ixxf), is_idx(xxf_lo,ixxf))
  end subroutine xxfidx2yxfidx

  elemental subroutine yxfidx2xxfidx (iky, iyxf, yxf_lo, xxf_lo, ikx, ixxf)
    implicit none
    integer, intent (in) :: iky, iyxf
    type (yxf_layout_type), intent (in) :: yxf_lo
    type (xxf_layout_type), intent (in) :: xxf_lo
    integer, intent (out) :: ikx, ixxf
    integer :: iky0

    iky0 = iky
    if (iky0 > xxf_lo%naky) then
       ikx = -999999
       ixxf = -999999
       return
    end if

    ikx = ikx_idx(yxf_lo,iyxf)
    ixxf = idx(xxf_lo, ig_idx(yxf_lo,iyxf), iv_idx(yxf_lo,iyxf), &
         iky0, imu_idx(yxf_lo,iyxf), is_idx(yxf_lo,iyxf))
  end subroutine yxfidx2xxfidx

  ! subroutine factors (n, j, div)
  !   integer, intent (in) :: n
  !   integer, intent (out) :: j
  !   integer, dimension (:), intent (out) :: div
  !   integer :: i, imax

  !   do i=2,n
  !      if (mod(n,i)==0) exit
  !   end do
  !   imax = n/i
  !   j=1
  !   do i=1,imax
  !      if (mod(n,i)==0) then
  !         div(j) = i
  !         j=j+1
  !      end if
  !   end do
  !   div(j) = n
  ! end subroutine factors

  subroutine finish_layouts
    
    implicit none
    
  end subroutine finish_layouts
  
end module stella_layouts
