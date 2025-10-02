!###############################################################################
!                                 ARRAYS FIELDS                                 
!###############################################################################
! This module stores the following fields:
!     - phi(kx,ky,z)
!     - apar(kx,ky,z)
!     - bpar(kx,ky,z)
! 
! So that all stella modules can easily accesss it.
! 
! The fields have dimensions: (naky, nakx, -nzgrid:nzgrid, ntubes)
! 
! All spatial information is local for every processor.
!###############################################################################
module arrays_fields

   use mpi
   use common_types, only: eigen_type

   implicit none
   
   ! Make the fields available to all modules
   public :: phi, phi_old
   public :: apar, apar_old
   public :: bpar, bpar_old
   public :: phi_shared
   public :: phi_corr_QN, apar_corr_QN
   public :: phi_proj, phi_proj_stage
   public :: phi_corr_GA, apar_corr_GA
   public :: phi_ext
   public :: phi_solve
   public :: phizf_solve
   
   private

   ! Electrostatic potential <phi> and electromagnetic fields <apar> and <bpar>
   ! (naky, nakx, -nzgrid:nzgrid, ntubes)
   complex, dimension(:, :, :, :), allocatable :: phi, phi_old
   complex, dimension(:, :, :, :), allocatable :: apar, apar_old
   complex, dimension(:, :, :, :), allocatable :: bpar, bpar_old

   ! DSO - the following is a band-aid for radially global simulations until we more fully incorporate shared memory
   ! (naky, nakx, -nzgrid:nzgrid, ntubes)
   complex, dimension(:, :, :, :), pointer :: phi_shared

   ! Radial corrections to phi and apar from quasineutrality/whatever controls apar
   ! (naky, nakx, -nzgrid:nzgrid, ntubes)
   complex, dimension(:, :, :, :), allocatable :: phi_corr_QN, apar_corr_QN

   ! Needed to implement time-delayed source when using projection method
   ! (nakx, -nzgrid:nzgrid, ntubes)
   complex, dimension(:, :, :), allocatable :: phi_proj, phi_proj_stage

   ! Radial corrections to phi and apar from gyroaveraging may result in tight space constraints
   ! (naky, nakx, -nzgrid:nzgrid, ntubes, -vmu-layout-)
   complex, dimension(:, :, :, :, :), allocatable :: phi_corr_GA, apar_corr_GA

   ! (nakx*nztot)
   complex, dimension(:), pointer :: phi_ext => null()

   ! TODO - dimensions?
   type(eigen_type), dimension(:, :), allocatable :: phi_solve
   type(eigen_type) :: phizf_solve

end module arrays_fields
