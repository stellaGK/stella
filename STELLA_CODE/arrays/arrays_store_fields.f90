module arrays_store_fields

   use mpi
   use stella_common_types, only: eigen_type

   implicit none

   complex, dimension(:, :, :, :), allocatable :: phi, phi_old
   complex, dimension(:, :, :, :), allocatable :: apar, apar_old
   complex, dimension(:, :, :, :), allocatable :: bpar, bpar_old
   ! (naky, nakx, -nzgrid:nzgrid, ntubes)

   ! DSO 0 the following is a band-aid for radially global simulations until
   ! we more fully incorporate shared memory
   complex, dimension(:, :, :, :), pointer :: phi_shared
   ! (naky, nakx, -nzgrid:nzgrid, ntubes)

   ! radial corrections to phi and apar from quasineutrality/whatever controls apar
   complex, dimension(:, :, :, :), allocatable :: phi_corr_QN, apar_corr_QN
   ! (naky, nakx, -nzgrid:nzgrid, ntubes)

   ! needed to implement time-delayed source when using projection method
   complex, dimension(:, :, :), allocatable :: phi_proj, phi_proj_stage
   ! (nakx, -nzgrid:nzgrid, ntubes)

   ! radial corrections to phi and apar from gyroaveraging
   ! may result in tight space constraints however
   complex, dimension(:, :, :, :, :), allocatable :: phi_corr_GA, apar_corr_GA
   ! (naky, nakx, -nzgrid:nzgrid, ntubes, -vmu-layout-)

   complex, dimension(:), pointer :: phi_ext => null()
   ! (nakx*nztot)

   type(eigen_type), dimension(:, :), allocatable :: phi_solve
   type(eigen_type) :: phizf_solve

end module arrays_store_fields
