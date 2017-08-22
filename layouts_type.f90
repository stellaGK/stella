module layouts_type
  ! This can be made by just replacing nakx by ntheta0
  !   from AstroGK's layouts_type.f90

  implicit none

  type :: kxkyz_layout_type
     sequence
     integer :: iproc
     integer :: ntgrid, nzed, naky, nakx, nvgrid, nvpa, nmu, nspec
     integer :: llim_world, ulim_world, llim_proc, ulim_proc, ulim_alloc, blocksize
  end type kxkyz_layout_type

  type :: kxyz_layout_type
     sequence
     integer :: iproc
     integer :: ntgrid, nzed, ny, naky, nakx, nvgrid, nvpa, nmu, nspec
     integer :: llim_world, ulim_world, llim_proc, ulim_proc, ulim_alloc, blocksize
  end type kxyz_layout_type

  type :: vmu_layout_type
     sequence
     logical :: xyz
     integer :: iproc
     integer :: ntgrid, nzed, ny, naky, nakx, nvgrid, nvpa, nmu, nspec
     integer :: llim_world, ulim_world, llim_proc, ulim_proc, ulim_alloc, blocksize
  end type vmu_layout_type

end module layouts_type
