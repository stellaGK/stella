module fields_arrays

  use common_types, only: response_matrix_type

  implicit none

  complex, dimension (:,:,:,:), allocatable :: phi, apar, phi_old
  ! (naky, nakx, -nzgrid:nzgrid, ntubes)

  type (response_matrix_type), dimension (:), allocatable :: response_matrix

end module fields_arrays
