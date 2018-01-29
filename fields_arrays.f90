module fields_arrays

  use common_types, only: response_matrix_type

  implicit none

  complex, dimension (:,:,:), allocatable :: phi, apar
  ! (naky, nakx, -nzgrid:nzgrid)

  complex, dimension (:,:), allocatable :: phi0_old
  ! (naky, nakx)

  type (response_matrix_type), dimension (:), allocatable :: response_matrix

end module fields_arrays
