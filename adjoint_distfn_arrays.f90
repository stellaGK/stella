module adjoint_distfn_arrays

   public :: g_omega, g_omega2
   public :: gsave
   public :: lam_save
   public :: g_unpert
   public :: g_store

   public :: source_adjoint

   complex, dimension(:, :, :, :, :), allocatable :: g_omega, g_omega2
   complex, dimension(:, :, :, :, :), allocatable :: gsave
   complex, dimension(:, :, :, :, :), allocatable :: lam_save

   complex, dimension(:, :, :, :, :), allocatable :: g_unpert

   complex, dimension(:, :, :, :, :), allocatable :: g_store

   complex, dimension(:, :, :, :, :), allocatable :: source_adjoint

end module adjoint_distfn_arrays

