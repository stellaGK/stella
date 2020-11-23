module fields_arrays

  use common_types, only: response_matrix_type

  implicit none

  complex, dimension (:,:,:,:), allocatable :: phi, apar, phi_old
  ! (naky, nakx, -nzgrid:nzgrid, ntubes)

  ! radial corrections to phi and apar from quasineutrality/whatever controls apar
  complex, dimension (:,:,:,:), allocatable :: phi_corr_QN, apar_corr_QN
  ! (naky, nakx, -nzgrid:nzgrid, ntubes)

  ! radial corrections to phi and apar from gyroaveraging
  ! may result in tight space constraints however
  complex, dimension (:,:,:,:,:), allocatable :: phi_corr_GA, apar_corr_GA
  complex, dimension (:,:,:,:,:), allocatable :: phi_corr_ztmax
  ! (naky, nakx, -nzgrid:nzgrid, ntubes, -vmu-layout-)

  type (response_matrix_type), dimension (:), allocatable :: response_matrix

  real, dimension (:), allocatable :: shift_state

end module fields_arrays
