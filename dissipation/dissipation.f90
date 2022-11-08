module dissipation

   implicit none

   public :: read_parameters
   public :: init_dissipation, finish_dissipation
   public :: init_collisions, collisions_initialized
   public :: advance_collisions_explicit, advance_collisions_implicit

   public :: include_collisions
   public :: hyper_dissipation
   public :: collisions_implicit
   public :: cfl_dt_mudiff, cfl_dt_vpadiff

   public :: time_collisions

   private

   logical :: collisions_initialized = .false.

   logical :: include_collisions
   logical :: collisions_implicit
   logical :: hyper_dissipation

   character(30) :: collision_model

   real :: cfl_dt_mudiff = huge(0.0), cfl_dt_vpadiff = huge(0.0)
   real, dimension(2, 2) :: time_collisions = 0.

contains

   subroutine init_dissipation

      use mp, only: proc0
      use hyper, only: init_hyper

      implicit none

      call read_parameters
      if (proc0) then
         write (*, '(A)') "############################################################"
         write (*, '(A)') "                         COLLISIONS"
         write (*, '(A)') "############################################################"
         if (include_collisions) then
            if (collision_model == "dougherty") then
               write (*, *)
               write (*, *) 'Coll. model:     Dougherty'
               if (collisions_implicit) then
                  write (*, *) 'Coll. algorithm: implicit'
               else
                  write (*, *) 'Coll. algorithm: explicit'
               end if
            else if (collision_model == "fokker-planck") then
               write (*, *) 'Coll. model:     Fokker-Planck'
               if (collisions_implicit) then
                  write (*, *) 'Coll. algorithm: implicit'
               else
                  write (*, *) 'Coll. algorithm: explicit'
               end if
            end if
            write (*, *)
         else
            write (*, *) 'Coll. model:     None'
            write (*, *)
         end if
      end if

      if (hyper_dissipation) then
         call init_hyper
      end if

   end subroutine init_dissipation

   subroutine read_parameters

      use file_utils, only: input_unit_exist
      use mp, only: proc0, broadcast
      use run_parameters, only: fully_explicit
      use coll_dougherty, only: read_parameters_dougherty
      use coll_fokkerplanck, only: read_parameters_fp
      use hyper, only: read_parameters_hyper

      implicit none

      namelist /dissipation/ include_collisions, collisions_implicit, collision_model, hyper_dissipation

      integer :: in_file
      logical :: dexist

      if (proc0) then
         include_collisions = .false.
         collisions_implicit = .true.
         collision_model = "dougherty"        ! dougherty or fokker-planck
         hyper_dissipation = .false.

         in_file = input_unit_exist("dissipation", dexist)
         if (dexist) read (unit=in_file, nml=dissipation)
      end if

      call broadcast(include_collisions)
      call broadcast(collisions_implicit)
      call broadcast(collision_model)
      call broadcast(hyper_dissipation)

      if (.not. include_collisions) collisions_implicit = .false.

      if (include_collisions) then
         if (collision_model == "dougherty") then
            call read_parameters_dougherty
         else if (collision_model == "fokker-planck") then
            call read_parameters_fp
         end if
      end if

      if (hyper_dissipation) then
         call read_parameters_hyper
         fully_explicit = .false.
      end if

   end subroutine read_parameters

   subroutine init_collisions

      use coll_dougherty, only: init_collisions_dougherty
      use coll_fokkerplanck, only: init_collisions_fp

      implicit none

      if (collisions_initialized) return
      collisions_initialized = .true.

      if (collision_model == "dougherty") then
         call init_collisions_dougherty(collisions_implicit, cfl_dt_vpadiff, cfl_dt_mudiff)
      else if (collision_model == "fokker-planck") then
         call init_collisions_fp(collisions_implicit, cfl_dt_vpadiff, cfl_dt_mudiff)
      end if

   end subroutine init_collisions

   subroutine finish_dissipation

      implicit none

      call finish_collisions

   end subroutine finish_dissipation

   subroutine finish_collisions

      use coll_dougherty, only: finish_collisions_dougherty
      use coll_fokkerplanck, only: finish_collisions_fp

      implicit none

      if (collision_model == "dougherty") then
         call finish_collisions_dougherty
      else if (collision_model == "fokker-planck") then
         call finish_collisions_fp
      end if

      collisions_initialized = .false.

   end subroutine finish_collisions

   subroutine advance_collisions_explicit(g, phi, gke_rhs)

      use mp, only: mp_abort
      use physics_flags, only: full_flux_surface
      use stella_layouts, only: vmu_lo
      use zgrid, only: nzgrid
      use coll_dougherty, only: advance_collisions_dougherty_explicit
      use coll_fokkerplanck, only: advance_collisions_fp_explicit

      implicit none

      complex, dimension(:, :, -nzgrid:, :, vmu_lo%llim_proc:), intent(in) :: g
      complex, dimension(:, :, -nzgrid:, :), intent(in) :: phi
      complex, dimension(:, :, -nzgrid:, :, vmu_lo%llim_proc:), intent(in out) :: gke_rhs

      if (full_flux_surface) then
         call mp_abort("collisions not currently supported for full_flux_surface=T.  Aborting.")
      end if

      if (collision_model == "dougherty") then
         call advance_collisions_dougherty_explicit(g, phi, gke_rhs, time_collisions)
      else if (collision_model == "fokker-planck") then
         call advance_collisions_fp_explicit(g, phi, gke_rhs, time_collisions)
      end if

   end subroutine advance_collisions_explicit

   subroutine advance_collisions_implicit(mirror_implicit, phi, apar, bpar, g)

      use mp, only: proc0
      use redistribute, only: gather, scatter
      use dist_redistribute, only: kxkyz2vmu
      use job_manage, only: time_message
      use zgrid, only: nzgrid
      use vpamu_grids, only: set_vpa_weights
      use stella_layouts, only: vmu_lo
      use dist_fn_arrays, only: gvmu
      use coll_dougherty, only: advance_collisions_dougherty_implicit
      use coll_fokkerplanck, only: advance_collisions_fp_implicit

      implicit none

      logical, intent(in) :: mirror_implicit
      complex, dimension(:, :, -nzgrid:, :), intent(in out) :: phi, apar, bpar
      complex, dimension(:, :, -nzgrid:, :, vmu_lo%llim_proc:), intent(in out) :: g

      logical :: conservative_wgts

      if (proc0) call time_message(.false., time_collisions(:, 1), ' collisions')

      !> switch the vpa integration weights to ensure correct integration by parts
      conservative_wgts = .true.
      call set_vpa_weights(conservative_wgts)

      if (proc0) call time_message(.false., time_collisions(:, 2), ' coll_redist')
      call scatter(kxkyz2vmu, g, gvmu)
      if (proc0) call time_message(.false., time_collisions(:, 2), ' coll_redist')

      if (collision_model == "dougherty") then
         call advance_collisions_dougherty_implicit(phi, apar, bpar)
      else if (collision_model == "fokker-planck") then
         call advance_collisions_fp_implicit(phi, apar, bpar)
      end if

      if (.not. mirror_implicit) then
         ! then take the results and remap again so ky,kx,z local.
         if (proc0) call time_message(.false., time_collisions(:, 2), ' coll_redist')
         call gather(kxkyz2vmu, gvmu, g)
         if (proc0) call time_message(.false., time_collisions(:, 2), ' coll_redist')
      end if

      conservative_wgts = .false.
      call set_vpa_weights(conservative_wgts)

      if (proc0) call time_message(.false., time_collisions(:, 1), ' collisions')

   end subroutine advance_collisions_implicit

end module dissipation
