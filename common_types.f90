module common_types

  implicit none
  
  type :: kxkyz_layout_type
     sequence
     integer :: iproc
     integer :: nzgrid, nzed, naky, nakx, nvgrid, nvpa, nmu, nspec
     integer :: llim_world, ulim_world, llim_proc, ulim_proc, ulim_alloc, blocksize
  end type kxkyz_layout_type
  
  type :: kxyz_layout_type
     sequence
     integer :: iproc
     integer :: nzgrid, nzed, ny, naky, nakx, nvgrid, nvpa, nmu, nspec
     integer :: llim_world, ulim_world, llim_proc, ulim_proc, ulim_alloc, blocksize
  end type kxyz_layout_type
  
  type :: vmu_layout_type
     sequence
     logical :: xyz
     integer :: iproc
     integer :: nzgrid, nzed, ny, naky, nakx, nvgrid, nvpa, nmu, nspec
     integer :: llim_world, ulim_world, llim_proc, ulim_proc, ulim_alloc, blocksize
  end type vmu_layout_type

  type :: flux_surface_type
     real :: rmaj
     real :: rgeo
     real :: kappa
     real :: kapprim
     real :: tri
     real :: triprim
     real :: rhoc
     real :: dr
     real :: shift
     real :: qinp
     real :: shat
     real :: betaprim
     real :: betadbprim
     real :: d2qdr2
     real :: d2psidr2
     real :: dpsitordrho
     real :: rhotor
     real :: drhotordrho
  end type flux_surface_type
  
  type spec_type
     real :: z
     real :: mass
     real :: dens, temp
     real :: tprim, fprim
     real :: vnewk, nustar, nu
     real :: stm, zstm, tz, smz, zt
     integer :: type
  end type spec_type
  
end module common_types
