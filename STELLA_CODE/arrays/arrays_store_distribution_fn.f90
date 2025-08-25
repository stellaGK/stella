module arrays_store_distribution_fn
   ! A container for the arrays that are used to store the distribution function arrays
   public :: gnew, gold, g_symm, g_scratch
   public :: g0, g1, g2, g3
   public :: g_krook, g_proj
   public :: gvmu, g_kymus

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

   ! gyro-averaged pdf used when advancing parallel streaming and mirror
   ! without operator splitting
   complex, dimension(:, :, :, :, :), allocatable :: g_kymus
   ! (nakx, -nzgrid:nzgrid, ntubes, vpa, -kymus-layout-)
   
   complex, dimension(:, :, :), allocatable :: gvmu
   ! (nvpa, nmu, -kxkyz-layout-)

end module arrays_store_distribution_fn
