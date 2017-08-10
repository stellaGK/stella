module layouts_type
  ! This can be made by just replacing nakx by ntheta0
  !   from AstroGK's layouts_type.f90

  implicit none

  type :: gvmu_layout_type
     sequence
     integer :: iproc
     integer :: ntgrid, ntheta, naky, ntheta0, nvgrid, nvpa, nmu, nspec
     integer :: llim_world, ulim_world, llim_proc, ulim_proc, ulim_alloc, blocksize
  end type gvmu_layout_type

  type :: gy_layout_type
     sequence
     integer :: iproc
     integer :: ntgrid, ntheta, naky, ntheta0, nvgrid, nvpa, nmu, nspec, ny
     integer :: llim_world, ulim_world, llim_proc, ulim_proc, ulim_alloc, blocksize
     integer :: it_ord, ig_ord, iv_ord, imu_ord, is_ord
     integer :: it_comp, ig_comp, iv_comp, imu_comp, is_comp
     integer, dimension (5) :: compound_count
  end type gy_layout_type

  type :: gxyz_layout_type
     sequence
     logical :: xyz
     integer :: iproc
     integer :: ntgrid, ntheta, naky, ntheta0, nvgrid, nvpa, nmu, nspec
     integer :: llim_world, ulim_world, llim_proc, ulim_proc, ulim_alloc, blocksize
  end type gxyz_layout_type

end module layouts_type
