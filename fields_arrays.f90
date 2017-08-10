module fields_arrays

  implicit none

  complex, dimension (:,:,:), allocatable :: phi
  complex, dimension (:,:,:), allocatable :: apar
!  complex, dimension (:), allocatable :: bparold, bparnew
  ! (gvmus_lo%llim_proc:gvmus_lo%ulim_alloc)

end module fields_arrays
