!###############################################################################
!                         ARRAYS DISTRIBUTION FUNCTION                          
!###############################################################################
! This module stores the guiding-center distribution function g(kx,ky,z,mu,vpa,s)
! so that all stella modules can easily accesss it.
!###############################################################################
module arrays_distribution_function

   implicit none

   ! Make the distribution functions available to all modules
   public :: gnew, gold, g_symm, g_scratch
   public :: g0, g1, g2, g3
   public :: g_krook, g_proj
   public :: gvmu, g_kymus
   
   private

   ! Distribution functions parallelised over (vpa, mu, s)
   ! (naky, nakx, -nzgrid:nzgrid, ntubes, -vmu-layout-)
   complex, dimension(:, :, :, :, :), allocatable :: gnew, gold, g_scratch
   complex, dimension(:, :, :, :, :), target, allocatable :: g_symm
   complex, dimension(:, :, :, :, :), allocatable :: g0, g1, g2, g3

   ! Needed to implement time-delayed source when using Krook operator
   ! (nakx, -nzgrid:nzgrid, ntubes, -vmu-layout-)
   complex, dimension(:, :, :, :), allocatable :: g_krook
   complex, dimension(:, :, :, :), allocatable :: g_proj

   ! Gyro-averaged pdf used when advancing parallel streaming and mirror without operator splitting
   ! Distribution function parallelised over (ky, mu, s)
   ! (nakx, -nzgrid:nzgrid, ntubes, vpa, -kymus-layout-)
   complex, dimension(:, :, :, :, :), allocatable :: g_kymus
   
   ! Distribution function parallelised over (kx,ky,z,s)
   ! (nvpa, nmu, -kxkyz-layout-)
   complex, dimension(:, :, :), allocatable :: gvmu

end module arrays_distribution_function
