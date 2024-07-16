module fields

   use common_types, only: eigen_type

   use mpi

   use common_types, only: coupled_alpha_type, gam0_ffs_type

   implicit none

   public :: init_fields, finish_fields
   public :: advance_fields

   public :: enforce_reality_field
   public :: rescale_fields
   public :: time_field_solve
   public :: fields_updated
   public :: get_dchidy, get_dchidx

   private

   logical :: fields_updated = .false.
   logical :: fields_initialized = .false.

   logical :: debug = .false.

   real, dimension(2) :: time_field_solve = 0.

   interface get_dchidy
      module procedure get_dchidy_4d
      module procedure get_dchidy_2d
   end interface get_dchidy

contains

   !> ######################################################################################
   !> Initialise fields arrays depending on type (i.e. Fluxtube, EM, FFS, Radially Global)
   !> ######################################################################################
   subroutine init_fields

      use mp, only: proc0
      use linear_solve, only: lu_decomposition
      use physics_flags, only: full_flux_surface, radial_variation
      use fields_fluxtube, only: init_fields_fluxtube
      use fields_ffs, only: init_fields_ffs
      use fields_radial_variation, only: allocate_arrays_radial_variation, init_radial_field_solve
      implicit none

      debug = debug .and. proc0
      if (fields_initialized) return
      fields_initialized = .true.

      !> Allocate arrays such as phi that are needed throughout the simulation
      call allocate_arrays

      if (full_flux_surface) then
         call init_fields_ffs
      else
         call init_fields_fluxtube
         if (radial_variation) then
            call allocate_arrays_radial_variation
            call init_radial_field_solve
         end if
      end if

      if (fields_initialized) return
      fields_initialized = .true.

   end subroutine init_fields

   !> Allocate arrays needed for solving fields for all versions of stella
   subroutine allocate_arrays

      use fields_arrays, only: phi, phi_old
      use fields_arrays, only: apar
      use fields_arrays, only: gamtot
      use zgrid, only: nzgrid, ntubes
      use stella_layouts, only: vmu_lo
      use kt_grids, only: naky, nakx

      implicit none

      if (.not. allocated(phi)) then
         allocate (phi(naky, nakx, -nzgrid:nzgrid, ntubes))
         phi = 0.
      end if

      if (.not. allocated(phi_old)) then
         allocate (phi_old(naky, nakx, -nzgrid:nzgrid, ntubes))
         phi_old = 0.
      end if

      !> TODO-GA: neeed to make this such that it is only for EM stella
      if (.not. allocated(apar)) then
         allocate (apar(naky, nakx, -nzgrid:nzgrid, ntubes))
         apar = 0.
      end if

      if (.not. allocated(gamtot)) then
         allocate (gamtot(naky, nakx, -nzgrid:nzgrid)); gamtot = 0.
      end if

   end subroutine allocate_arrays

   subroutine enforce_reality_field(fin)

!DSO> While most of the modes in the box have reality built in (as we
!     throw out half the kx-ky plane, modes with ky=0 do not have
!     this enforcement built in. In theory this should not be a problem
!     as these modes should be stable, but I made this function (and
!     its relative in the dist file) just in case

      use kt_grids, only: nakx
      use zgrid, only: nzgrid

      implicit none

      complex, dimension(:, :, -nzgrid:, :), intent(inout) :: fin

      integer ikx

      fin(1, 1, :, :) = real(fin(1, 1, :, :))
      do ikx = 2, nakx / 2 + 1
         fin(1, ikx, :, :) = 0.5 * (fin(1, ikx, :, :) + conjg(fin(1, nakx - ikx + 2, :, :)))
         fin(1, nakx - ikx + 2, :, :) = conjg(fin(1, ikx, :, :))
      end do

   end subroutine enforce_reality_field

   !> ######################################################################################
   !> Advance fields arrays during main code
   !> ######################################################################################
   subroutine advance_fields(g, phi, apar, dist)

      use mp, only: proc0
      use stella_layouts, only: vmu_lo
      use job_manage, only: time_message
      use redistribute, only: scatter
      use dist_fn_arrays, only: gvmu
      use zgrid, only: nzgrid
      use dist_redistribute, only: kxkyz2vmu
      use run_parameters, only: fields_kxkyz
      use physics_flags, only: full_flux_surface

      use fields_fluxtube, only: advance_fields_fluxtube
      use fields_arrays, only: time_field_solve

      implicit none

      complex, dimension(:, :, -nzgrid:, :, vmu_lo%llim_proc:), intent(in) :: g
      complex, dimension(:, :, -nzgrid:, :), intent(out) :: phi, apar
      character(*), intent(in) :: dist

      if (fields_updated) return

      !> time the communications + field solve
      if (proc0) call time_message(.false., time_field_solve(:, 1), ' fields')
      !> fields_kxkyz = F is the default

      !> Get phi for fluxtube
      call advance_fields_fluxtube(g, phi, apar, dist)

      !> set a flag to indicate that the fields have been updated
      !> this helps avoid unnecessary field solves
      fields_updated = .true.
      !> time the communications + field solve
      if (proc0) call time_message(.false., time_field_solve(:, 1), ' fields')

   end subroutine advance_fields

   !> ######################################################################################
   !> Routines that use fields
   !> ######################################################################################
   !> rescale fields, including the distribution function
   subroutine rescale_fields(target_amplitude)

      use mp, only: scope, subprocs, crossdomprocs, sum_allreduce
      use fields_arrays, only: phi, apar
      use dist_fn_arrays, only: gnew, gvmu
      use volume_averages, only: volume_average
      use job_manage, only: njobs
      use file_utils, only: runtype_option_switch, runtype_multibox

      implicit none

      real, intent(in) :: target_amplitude
      real :: phi2, rescale

      call volume_average(phi, phi2)

      if (runtype_option_switch == runtype_multibox) then
         call scope(crossdomprocs)
         call sum_allreduce(phi2)
         call scope(subprocs)
         phi2 = phi2 / njobs
      end if

      rescale = target_amplitude / sqrt(phi2)

      phi = rescale * phi
      apar = rescale * apar
      gnew = rescale * gnew
      gvmu = rescale * gvmu

   end subroutine rescale_fields

   !> compute d<chi>/dy in (ky,kx,z,tube) space
   subroutine get_dchidy_4d(phi, apar, dchidy)

      use constants, only: zi
      use gyro_averages, only: gyro_average
      use stella_layouts, only: vmu_lo
      use stella_layouts, only: is_idx, iv_idx
      use run_parameters, only: fphi, fapar
      use species, only: spec
      use zgrid, only: nzgrid, ntubes
      use vpamu_grids, only: vpa
      use kt_grids, only: nakx, aky, naky

      implicit none

      complex, dimension(:, :, -nzgrid:, :), intent(in) :: phi, apar
      complex, dimension(:, :, -nzgrid:, :, vmu_lo%llim_proc:), intent(out) :: dchidy

      integer :: ivmu, iv, is, iz, it, ikx
      complex, dimension(:, :, :, :), allocatable :: field

      allocate (field(naky, nakx, -nzgrid:nzgrid, ntubes))

      do ivmu = vmu_lo%llim_proc, vmu_lo%ulim_proc
         is = is_idx(vmu_lo, ivmu)
         iv = iv_idx(vmu_lo, ivmu)
         do it = 1, ntubes
            do iz = -nzgrid, nzgrid
               do ikx = 1, nakx
                  field(:, ikx, iz, it) = zi * aky(:) * (fphi * phi(:, ikx, iz, it) - fapar * vpa(iv) * spec(is)%stm * apar(:, ikx, iz, it))
               end do
            end do
         end do
         call gyro_average(field, ivmu, dchidy(:, :, :, :, ivmu))
      end do

      deallocate (field)

   end subroutine get_dchidy_4d

   !> compute d<chi>/dy in (ky,kx) space
   subroutine get_dchidy_2d(iz, ivmu, phi, apar, dchidy)

      use constants, only: zi
      use gyro_averages, only: gyro_average
      use stella_layouts, only: vmu_lo
      use stella_layouts, only: is_idx, iv_idx
      use run_parameters, only: fphi, fapar
      use species, only: spec
      use vpamu_grids, only: vpa
      use kt_grids, only: nakx, aky, naky

      implicit none

      integer, intent(in) :: ivmu, iz
      complex, dimension(:, :), intent(in) :: phi, apar
      complex, dimension(:, :), intent(out) :: dchidy

      integer :: iv, is
      complex, dimension(:, :), allocatable :: field

      allocate (field(naky, nakx))

      is = is_idx(vmu_lo, ivmu)
      iv = iv_idx(vmu_lo, ivmu)
      field = zi * spread(aky, 2, nakx) &
              * (fphi * phi - fapar * vpa(iv) * spec(is)%stm * apar)
      call gyro_average(field, iz, ivmu, dchidy)

      deallocate (field)

   end subroutine get_dchidy_2d

   !> compute d<chi>/dx in (ky,kx) space
   subroutine get_dchidx(iz, ivmu, phi, apar, dchidx)

      use constants, only: zi
      use gyro_averages, only: gyro_average
      use stella_layouts, only: vmu_lo
      use stella_layouts, only: is_idx, iv_idx
      use run_parameters, only: fphi, fapar
      use species, only: spec
      use vpamu_grids, only: vpa
      use kt_grids, only: akx, naky, nakx

      implicit none

      integer, intent(in) :: ivmu, iz
      complex, dimension(:, :), intent(in) :: phi, apar
      complex, dimension(:, :), intent(out) :: dchidx

      integer :: iv, is
      complex, dimension(:, :), allocatable :: field

      allocate (field(naky, nakx))

      is = is_idx(vmu_lo, ivmu)
      iv = iv_idx(vmu_lo, ivmu)
      field = zi * spread(akx, 1, naky) &
              * (fphi * phi - fapar * vpa(iv) * spec(is)%stm * apar)
      call gyro_average(field, iz, ivmu, dchidx)

      deallocate (field)

   end subroutine get_dchidx

   !> ######################################################################################
   !> Deallocate arrays for fields
   !> ######################################################################################

   subroutine finish_fields

      use fields_arrays, only: phi, phi_old
      use fields_arrays, only: gamtot, gamtot3
      use fields_arrays, only: c_mat, theta
      use physics_flags, only: full_flux_surface, radial_variation
      use fields_ffs, only: finish_fields_ffs
      use fields_radial_variation, only: finish_radial_fields
      !> TODO-GA: move apar stuff to EM fields
      use fields_arrays, only: apar
      implicit none

      if (allocated(phi)) deallocate (phi)
      if (allocated(phi_old)) deallocate (phi_old)
      if (allocated(gamtot)) deallocate (gamtot)
      if (allocated(gamtot3)) deallocate (gamtot3)

      if (allocated(c_mat)) deallocate (c_mat)
      if (allocated(theta)) deallocate (theta)

      !> TODO-GA: REMOVE
      if (allocated(apar)) deallocate (apar)
      if (full_flux_surface) call finish_fields_ffs
      if (radial_variation) call finish_radial_fields

      fields_initialized = .false.

   end subroutine finish_fields

end module fields
