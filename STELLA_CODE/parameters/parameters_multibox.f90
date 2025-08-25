module parameters_multibox

   implicit none

   public :: read_multibox_parameters

   public :: ky_solve_real, ky_solve_radial
   public :: include_pressure_variation 
   public :: include_geometric_variation
   public :: smooth_zf
   public :: lr_debug_switch, lr_debug_option_default, & 
            lr_debug_option_L, lr_debug_option_R
   public :: krook_option_switch, krook_option_default, &
            krook_option_flat, krook_option_linear, &
            krook_option_exp, krook_option_exp_rev
   public :: mb_zf_option_switch, mb_zf_option_default, &
            mb_zf_option_skip_ky0, mb_zf_option_zero_ky0, &
            mb_zf_option_zero_fsa
   public :: rk_step
   public :: nu_krook_mb
   public :: mb_debug_step
   public :: krook_exponent, krook_efold
   public :: phi_bound, phi_pow
   public :: boundary_size, krook_size
   public :: include_multibox_krook
   public :: use_dirichlet_bc

   private

   logical :: ky_solve_real
   integer :: ky_solve_radial
   logical :: include_pressure_variation 
   logical :: include_geometric_variation
   logical :: smooth_zf

   integer :: lr_debug_switch
   integer, parameter:: lr_debug_option_default = 0, &
                     lr_debug_option_L = 1, &
                     lr_debug_option_R = 2

   integer :: krook_option_switch
   integer, parameter:: krook_option_default = 2, &
                           krook_option_flat = 0, &
                           krook_option_linear = 1, &
                           krook_option_exp = 2, &
                           krook_option_exp_rev = 3

   integer:: mb_zf_option_switch
   integer, parameter :: mb_zf_option_default = 0, &
                           mb_zf_option_skip_ky0 = 1, &
                           mb_zf_option_zero_ky0 = 2, &
                           mb_zf_option_zero_fsa = 3

   logical :: rk_step
   real :: nu_krook_mb
   integer :: mb_debug_step
   real :: krook_exponent, krook_efold
   real :: phi_bound, phi_pow
   logical :: use_dirichlet_BC
   integer :: boundary_size, krook_size
   logical :: include_multibox_krook
   
contains

   subroutine read_multibox_parameters

      use mp, only: broadcast
      use namelist_radial_variation, only: read_namelist_radial_variation

      implicit none

      call read_namelist_radial_variation(ky_solve_real, ky_solve_radial, &
               include_pressure_variation, include_geometric_variation, &
               smooth_zf, lr_debug_switch, krook_option_switch, mb_zf_option_switch, &
               rk_step, nu_krook_mb, mb_debug_step, &
               krook_exponent, krook_efold, phi_bound, phi_pow, &
               use_dirichlet_bc, boundary_size, krook_size)

      call broadcast(boundary_size)
      call broadcast(krook_size)
      call broadcast(nu_krook_mb)
      call broadcast(smooth_zf)
      call broadcast(mb_zf_option_switch)
      call broadcast(krook_option_switch)
      call broadcast(krook_exponent)
      call broadcast(krook_efold)
      call broadcast(lr_debug_switch)
      call broadcast(rk_step)
      call broadcast(mb_debug_step)
      call broadcast(phi_bound)
      call broadcast(phi_pow)
      call broadcast(use_dirichlet_BC)
      call broadcast(ky_solve_radial)
      call broadcast(ky_solve_real)
      call broadcast(include_pressure_variation)
      call broadcast(include_geometric_variation)
      call broadcast(use_dirichlet_bc)

   end subroutine read_multibox_parameters

end module parameters_multibox