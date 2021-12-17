module fields_arrays

   use common_types, only: response_matrix_type, eigen_type

   implicit none

   complex, dimension(:, :, :, :), allocatable :: phi, apar, bpar, phi_old
   ! (naky, nakx, -nzgrid:nzgrid, ntubes)

   ! radial corrections to phi and apar from quasineutrality/whatever controls apar
   complex, dimension(:, :, :, :), allocatable :: phi_corr_QN, apar_corr_QN
   ! (naky, nakx, -nzgrid:nzgrid, ntubes)

   ! radial corrections to phi and apar from gyroaveraging
   ! may result in tight space constraints however
   complex, dimension(:, :, :, :, :), allocatable :: phi_corr_GA, apar_corr_GA
   complex, dimension(:, :, :, :, :), allocatable :: phi_corr_ztmax
   ! (naky, nakx, -nzgrid:nzgrid, ntubes, -vmu-layout-)

   type(response_matrix_type), dimension(:), allocatable :: response_matrix

   real, dimension(:), allocatable :: shift_state

   real, dimension(:, :, :), allocatable :: gamtot, dgamtotdr, gamtot13, gamtot31, gamtot33, apar_denom !, gamone

   real :: gamtot_h, gamtot3_h, efac, efacp

   type(eigen_type), dimension(:, :), allocatable :: phi_solve
   type(eigen_type) :: phizf_solve

end module fields_arrays
