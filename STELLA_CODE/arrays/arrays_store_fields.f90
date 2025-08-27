!###############################################################################
!                                                                               
!###############################################################################
! This module ...
!###############################################################################
module arrays_store_fields

   use mpi
   use stella_common_types, only: eigen_type

   implicit none

   ! (naky, nakx, -nzgrid:nzgrid, ntubes)
   complex, dimension(:, :, :, :), allocatable :: phi, phi_old
   complex, dimension(:, :, :, :), allocatable :: apar, apar_old
   complex, dimension(:, :, :, :), allocatable :: bpar, bpar_old

   ! DSO 0 the following is a band-aid for radially global simulations until
   ! we more fully incorporate shared memory
   ! (naky, nakx, -nzgrid:nzgrid, ntubes)
   complex, dimension(:, :, :, :), pointer :: phi_shared

   ! Radial corrections to phi and apar from quasineutrality/whatever controls apar
   ! (naky, nakx, -nzgrid:nzgrid, ntubes)
   complex, dimension(:, :, :, :), allocatable :: phi_corr_QN, apar_corr_QN

   ! Needed to implement time-delayed source when using projection method
   ! (nakx, -nzgrid:nzgrid, ntubes)
   complex, dimension(:, :, :), allocatable :: phi_proj, phi_proj_stage

   ! Radial corrections to phi and apar from gyroaveraging
   ! may result in tight space constraints however
   ! (naky, nakx, -nzgrid:nzgrid, ntubes, -vmu-layout-)
   complex, dimension(:, :, :, :, :), allocatable :: phi_corr_GA, apar_corr_GA

   ! (nakx*nztot)
   complex, dimension(:), pointer :: phi_ext => null()

   ! TODO - dimensions?
   type(eigen_type), dimension(:, :), allocatable :: phi_solve
   type(eigen_type) :: phizf_solve

end module arrays_store_fields
