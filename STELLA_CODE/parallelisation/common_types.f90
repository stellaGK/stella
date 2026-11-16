!###############################################################################
!                                                                               
!###############################################################################
! This module ...
!###############################################################################
module common_types

   implicit none

   type :: kxkyz_layout_type
      sequence
      integer :: iproc
      integer :: nzgrid, nzed, ntubes, naky, nakx, nvgrid, nvpa, nmu, nspec
      integer :: llim_world, ulim_world, llim_proc, ulim_proc, ulim_alloc, blocksize
   end type kxkyz_layout_type

   type :: kxyz_layout_type
      sequence
      integer :: iproc
      integer :: nzgrid, nzed, ntubes, ny, naky, nakx, ikx_max, nvgrid, nvpa, nmu, nspec
      integer :: llim_world, ulim_world, llim_proc, ulim_proc, ulim_alloc, blocksize
   end type kxyz_layout_type

   type :: xyz_layout_type
      sequence
      integer :: iproc
      integer :: nzgrid, nzed, ntubes, ny, naky, nx, nakx, nvgrid, nvpa, nmu, nspec
      integer :: llim_world, ulim_world, llim_proc, ulim_proc, ulim_alloc, blocksize
   end type xyz_layout_type

   type :: vmu_layout_type
      sequence
      logical :: xyz
      integer :: iproc
      integer :: nzgrid, nzed, ntubes, nalpha, ny, naky, nx, nakx, nvgrid, nvpa, nmu, nspec
      integer :: llim_world, ulim_world, llim_proc, ulim_proc, ulim_alloc, blocksize
   end type vmu_layout_type

   type :: kymus_layout_type
      sequence
      integer :: iproc
      integer :: nzgrid, nzed, ntubes, naky, nakx, nvgrid, nvpa, nmu, nspec
      integer :: llim_world, ulim_world, llim_proc, ulim_proc, ulim_alloc, blocksize
   end type kymus_layout_type
   
   ! Not all quantities are calculated for every geometry (e.g. the Miller parameters and the radial
   ! derivatives for VMEC), so initialise them to zero to avoid using uninitialised values
   type :: flux_surface_type
      real :: rmaj = 0.0
      real :: rgeo = 0.0
      real :: kappa = 0.0
      real :: kapprim = 0.0
      real :: tri = 0.0
      real :: triprim = 0.0
      real :: rhoc = 0.0
      real :: dr = 0.0
      real :: shift = 0.0
      real :: qinp = 0.0
      real :: shat = 0.0
      real :: betaprim = 0.0
      real :: betadbprim = 0.0
      real :: d2qdr2 = 0.0
      real :: d2psidr2 = 0.0
      real :: dpsitordrho = 0.0
      real :: d2psitordrho2 = 0.0
      real :: rhotor = 0.0
      real :: drhotordrho = 0.0
      real :: psitor_lcfs = 0.0
      real :: zed0_fac = 0.0
      
      ! The next few variables are for multibox simulations
      real :: rhoc_psi0 = 0.0
      real :: qinp_psi0 = 0.0
      real :: shat_psi0 = 0.0
   end type flux_surface_type

   type spec_type
      integer :: nspec
      real :: z
      real :: mass
      real :: dens, temp
      real :: tprim, fprim
      real :: vnew_ref
      real :: stm, zstm, tz, smz, zt
      real :: d2ndr2, d2Tdr2
      real :: bess_fac ! 0 for argument of Bessel Functions equal to 0
      ! pre-2003 Fortran does not support
      ! allocatable arrays within derived types
      ! so set size large enough that it should be a problem
      ! should be nspec large
      real, dimension(10) :: vnew
      integer :: type

      ! The next few variables are for multibox simulations
      real :: dens_psi0, temp_psi0
      real :: stm_psi0, zstm_psi0, tz_psi0, smz_psi0, zt_psi0
   end type spec_type

   type :: eigen_type
      complex, dimension(:, :), pointer :: zloc => null()
      integer, dimension(:), pointer :: idx => null()
   end type eigen_type

   type :: response_matrix_type
      type(eigen_type), dimension(:), pointer :: eigen => null()
   end type response_matrix_type

   type coupled_alpha_type
      integer :: max_idx
      complex, dimension(:), pointer :: fourier => null()
   end type coupled_alpha_type

   type gam0_ffs_type
      integer, dimension(:), pointer :: pivot_index => null()
      complex, dimension(:, :), pointer :: matrix => null()
   end type gam0_ffs_type

end module common_types
