module adjoint_field_arrays

   implicit none

   complex, dimension(:, :, :, :), allocatable :: phi_save, chi_save

   complex, dimension(:, :), allocatable :: denominator
   complex, dimension(:, :), allocatable :: derivative

   complex, dimension(:, :), allocatable :: omega_g, omega

   complex, dimension(:, :, :, :), allocatable :: q_unpert
   complex, dimension(:, :, :, :), allocatable :: q_store

end module adjoint_field_arrays
