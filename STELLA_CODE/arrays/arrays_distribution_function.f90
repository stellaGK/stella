!###############################################################################
!                         ARRAYS DISTRIBUTION FUNCTION                          
!###############################################################################
! This module stores the guiding-center distribution function g(kx,ky,z,mu,vpa,s)
! so that all stella modules can easily accesss it.
! There are two possible layouts for the distribution function:
! 1) (naky, nakx, -nzgrid:nzgrid, ntubes, -vmu-layout-)
!    This is the more common layout that is used in the time advance,
!    and has all the spacial information local, whilst we parallelise
!    over the velocity space and species.
!    (used for the nonlinear term, drift advances, 
!    parallel streaming etc.)
! 2) (nvpa, nmu, -kxkyz-layout-)
!    This layout is used for anything that requires the velocity distibution to 
!    be local, such as in the mirror advance, or for collisions and the
!    calculation of moments.
! 3) (nakx, -nzgrid:nzgrid, ntubes, vpa, -kymus-layout-)
!    This layout is used for the parallel streaming and mirror advance when
!    not using operator splitting, as we need to have both the spatial and
!    vpa information local.
!###############################################################################
module arrays_distribution_function

   implicit none

   ! Make the distribution functions available to all modules
   public :: gnew, gold, g_symm, phi_gyro
   public :: g0, g1, g2, g3
   public :: g_krook, g_proj
   public :: gvmu, g_kymus
   
   private

   ! Distribution functions parallelised over (vpa, mu, s)
   ! (naky, nakx, -nzgrid:nzgrid, ntubes, -vmu-layout-)
   complex, dimension(:, :, :, :, :), allocatable :: gnew, gold, phi_gyro
   complex, dimension(:, :, :, :, :), target, allocatable :: g_symm
   complex, dimension(:, :, :, :, :), allocatable :: g0, g1, g2, g3

   ! Needed to implement time-delayed source when using Krook operator
   ! (nakx, -nzgrid:nzgrid, ntubes, -vmu-layout-)
   complex, dimension(:, :, :, :), allocatable :: g_krook
   complex, dimension(:, :, :, :), allocatable :: g_proj

   ! Gyro-averaged pdf used when advancing parallel streaming and mirror without 
   ! operator splitting. Distribution function is parallelised over (ky, mu, s)
   ! (nakx, -nzgrid:nzgrid, ntubes, vpa, -kymus-layout-)
   complex, dimension(:, :, :, :, :), allocatable :: g_kymus
   
   ! Distribution function parallelised over (kx,ky,z,s)
   ! (nvpa, nmu, -kxkyz-layout-)
   complex, dimension(:, :, :), allocatable :: gvmu

end module arrays_distribution_function
