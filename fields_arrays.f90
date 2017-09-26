module fields_arrays

  implicit none

  complex, dimension (:,:,:), allocatable :: phi, phi_old
  complex, dimension (:,:,:), allocatable :: apar, apar_old
  ! (naky, nakx, -nzgrid:nzgrid)

end module fields_arrays
