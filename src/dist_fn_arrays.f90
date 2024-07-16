!> A container for the arrays that are used to store the distribution function among other things.
!!  These need to be accessible at a lower dependency level than the dist_fn module itself.
!! These arrays are allocated in the function dist_fn::allocate_arrays.

module dist_fn_arrays

   public :: gnew, gold, g_symm, g_scratch
   public :: g0, g1, g2, g3
   public :: g_krook, g_proj
   public :: gvmu
   public :: kperp2, dkperp2dr
   public :: wstar, wstarp
   public :: wdriftx_g, wdrifty_g
   public :: wdriftx_phi, wdrifty_phi
   public :: wdriftx_bpar, wdrifty_bpar
   public :: wdriftpx_g, wdriftpy_g
   public :: wdriftpx_phi, wdriftpy_phi

   ! dist fn
   complex, dimension(:, :, :, :, :), allocatable :: gnew, gold, g_scratch
   ! (naky, nakx, -nzgrid:nzgrid, ntubes, -vmu-layout-)

   complex, dimension(:, :, :, :, :), target, allocatable :: g_symm
   ! (naky, nakx, -nzgrid:nzgrid, ntubes, -vmu-layout-)

   complex, dimension(:, :, :, :, :), allocatable :: g0, g1, g2, g3
   ! (naky, nakx, -nzgrid:nzgrid, ntubes, -vmu-layout-)

   ! needed to implement time-delayed source when using Krook operator
   complex, dimension(:, :, :, :), allocatable :: g_krook
   ! (nakx, -nzgrid:nzgrid, ntubes, -vmu-layout-)

   ! needed to implement time-delayed source when using projection method
   complex, dimension(:, :, :, :), allocatable :: g_proj
   ! (nakx, -nzgrid:nzgrid, ntubes, -vmu-layout-)

   complex, dimension(:, :, :), allocatable :: gvmu
   ! (nvpa, nmu, -kxkyz-layout-)

   real, dimension(:, :, :), allocatable :: wstar, wstarp
   ! (nalpha, -nzgrid:nzgrid, -vmu-layout-)

   real, dimension(:, :, :), allocatable :: wdriftx_g, wdrifty_g
   real, dimension(:, :, :), allocatable :: wdriftx_phi, wdrifty_phi
   real, dimension(:, :, :), allocatable :: wdriftx_bpar, wdrifty_bpar

   real, dimension(:, :, :), allocatable :: wdriftpx_g, wdriftpy_g
   real, dimension(:, :, :), allocatable :: wdriftpx_phi, wdriftpy_phi
   ! (nalpha, -nzgrid:nzgrid, -vmu-layout-)

   !> dkperp2dr will contain the radial variation of kperp2
   real, dimension(:, :, :, :), allocatable :: kperp2, dkperp2dr
   ! (naky, nakx, nalpha, -nzgrid:nzgrid)
   ! note: dkperp2dr is divided by kperp2

end module dist_fn_arrays
