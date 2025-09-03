module parameters_multibox

   ! Read the parameters for <lr_debug_switch> from namelist_radial_variation.f90
   use namelist_radial_variation, only: lr_debug_option_default
   use namelist_radial_variation, only: lr_debug_option_L
   use namelist_radial_variation, only: lr_debug_option_R
   
   ! Read the parameters for <krook_option_switch> from namelist_radial_variation.f90
   use namelist_radial_variation, only: krook_option_default
   use namelist_radial_variation, only: krook_option_flat
   use namelist_radial_variation, only: krook_option_linear
   use namelist_radial_variation, only: krook_option_exp
   use namelist_radial_variation, only: krook_option_exp_rev
   
   ! Read the parameters for <mb_zf_option_switch> from namelist_radial_variation.f90
   use namelist_radial_variation, only: mb_zf_option_default
   use namelist_radial_variation, only: mb_zf_option_skip_ky0
   use namelist_radial_variation, only: mb_zf_option_zero_ky0
   use namelist_radial_variation, only: mb_zf_option_zero_fsa

   implicit none

   ! Although the parameters are available through namelist_species, 
   ! make them available through grids_species as well
   public :: lr_debug_switch, lr_debug_option_default
   public :: lr_debug_option_L, lr_debug_option_R
   public :: krook_option_switch, krook_option_default
   public :: krook_option_flat, krook_option_linear
   public :: krook_option_exp, krook_option_exp_rev
   public :: mb_zf_option_switch, mb_zf_option_default
   public :: mb_zf_option_skip_ky0, mb_zf_option_zero_ky0
   public :: mb_zf_option_zero_fsa

   ! Make routines available to other modules
   public :: read_parameters_multibox
   public :: ky_solve_real, ky_solve_radial
   public :: include_pressure_variation 
   public :: include_geometric_variation
   public :: smooth_zf
   public :: rk_step
   public :: nu_krook_mb
   public :: mb_debug_step
   public :: krook_exponent, krook_efold
   public :: phi_bound, phi_pow
   public :: boundary_size, krook_size
   public :: include_multibox_krook
   public :: use_dirichlet_bc

   private

   ! Text option switches
   integer :: lr_debug_switch
   integer :: krook_option_switch
   integer:: mb_zf_option_switch
   
   logical :: ky_solve_real
   integer :: ky_solve_radial
   logical :: include_pressure_variation 
   logical :: include_geometric_variation
   logical :: smooth_zf
   logical :: rk_step
   real :: nu_krook_mb
   integer :: mb_debug_step
   real :: krook_exponent, krook_efold
   real :: phi_bound, phi_pow
   logical :: use_dirichlet_BC
   integer :: boundary_size, krook_size
   logical :: include_multibox_krook
   
contains

   !****************************************************************************
   !                                      Title
   !****************************************************************************
   subroutine read_parameters_multibox

      ! Parallelisation
      use mp, only: broadcast
      
      ! Read namelists from input file
      use namelist_radial_variation, only: read_namelist_multibox

      implicit none

      !-------------------------------------------------------------------------

      ! Read the "multibox_parameters" namelist in the input file
      call read_namelist_multibox(ky_solve_real, ky_solve_radial, &
         include_pressure_variation, include_geometric_variation, &
         smooth_zf, lr_debug_switch, krook_option_switch, mb_zf_option_switch, &
         rk_step, nu_krook_mb, mb_debug_step, &
         krook_exponent, krook_efold, phi_bound, phi_pow, &
         use_dirichlet_bc, boundary_size, krook_size)

      ! Broadcast the input parameters
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

   end subroutine read_parameters_multibox

end module parameters_multibox
