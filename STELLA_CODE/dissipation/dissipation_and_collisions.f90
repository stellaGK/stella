module dissipation

   implicit none

   ! Make routines available to other modules
   public :: read_parameters_dissipation_and_collisions
   public :: init_dissipation, finish_dissipation
   public :: init_collisions, initialised_collisions
   public :: advance_collisions_explicit, advance_collisions_implicit
   
   ! Make parmeters available
   public :: include_collisions, hyper_dissipation, ecoll_zeff
   public :: collisions_implicit, cfl_dt_mudiff, cfl_dt_vpadiff
   public :: time_collisions

   private

   ! Flags
   logical :: include_collisions
   logical :: collisions_implicit
   logical :: hyper_dissipation
   logical :: ecoll_zeff
   
   ! Collision model is "dougherty" or "fokker-planck"
   character(30) :: collision_model

   ! Parameters
   real :: cfl_dt_mudiff = huge(0.0), cfl_dt_vpadiff = huge(0.0)
   real, dimension(2, 2) :: time_collisions = 0.
   
   ! Only initialise once
   logical :: initialised_read_dissipation_and_collisions = .false.
   logical :: initialised_dissipation = .false.
   logical :: initialised_collisions = .false.

contains

!###############################################################################
!################################ READ NAMELIST ################################
!###############################################################################

   !****************************************************************************
   subroutine read_parameters_dissipation_and_collisions

      use mp, only: proc0, broadcast
      use file_utils, only: input_unit_exist
      use debug_flags, only: print_extra_info_to_terminal
      
      ! Read the <dissipation> namelist in the input file
      use namelist_dissipation, only: read_namelist_dissipation_and_collisions_options
      
      ! Read other input parameters related to specific collision models
      use dissipation_coll_dougherty, only: read_parameters_dougherty
      use dissipation_coll_fokkerplanck, only: read_parameters_fp
      use dissipation_hyper, only: read_parameters_hyper
      use parameters_numerical, only: fully_explicit

      implicit none
      
      !-------------------------------------------------------------------------

      ! Only initialise once
      if (initialised_read_dissipation_and_collisions) return
      initialised_read_dissipation_and_collisions = .true.
      
      ! Read <dissipation> namelist in the input file
      if (proc0) call read_namelist_dissipation_and_collisions_options(include_collisions, &
         collisions_implicit, collision_model, hyper_dissipation, ecoll_zeff)

      ! Broadcast to all processors
      call broadcast(include_collisions)
      call broadcast(collisions_implicit)
      call broadcast(collision_model)
      call broadcast(hyper_dissipation)
      call broadcast(ecoll_zeff)

      ! Read input parameters for the dougherty or fokker-planck collision model
      if (include_collisions) then
         if (collision_model == "dougherty") then
            call read_parameters_dougherty
         else if (collision_model == "fokker-planck") then
            call read_parameters_fp
         end if
      end if

      ! Read input parameters for the hyper-dissipation
      if (hyper_dissipation) then
         call read_parameters_hyper
         fully_explicit = .false.
      end if

      ! Print collision information to terminal
      if (proc0 .and. print_extra_info_to_terminal) then
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

   end subroutine read_parameters_dissipation_and_collisions
   
   
!###############################################################################
!#################### INITIALISE DISSIPATION AND COLLISIONS #################### 
!###############################################################################

   !--------------------------- Initialise dissipation -------------------------
   subroutine init_dissipation
   
      use dissipation_hyper, only: init_hyper

      implicit none
      
      !-------------------------------------------------------------------------

      ! Only initialise once
      if (initialised_dissipation) return
      initialised_dissipation = .true.

      ! Initialise hyper dissipation
      if (hyper_dissipation) then
         call init_hyper
      end if
   
   end subroutine init_dissipation

   !--------------------------- Initialise collisions --------------------------
   subroutine init_collisions

      use dissipation_coll_dougherty, only: init_collisions_dougherty
      use dissipation_coll_fokkerplanck, only: init_collisions_fp

      implicit none
      
      !-------------------------------------------------------------------------

      ! Only initialise once
      if (initialised_collisions) return
      initialised_collisions = .true.

      if (collision_model == "dougherty") then
         call init_collisions_dougherty(collisions_implicit, cfl_dt_vpadiff, cfl_dt_mudiff)
      else if (collision_model == "fokker-planck") then
         call init_collisions_fp(collisions_implicit, cfl_dt_vpadiff, cfl_dt_mudiff)
      end if

   end subroutine init_collisions

   !****************************************************************************
   subroutine finish_dissipation

      implicit none

      call finish_collisions

   end subroutine finish_dissipation

   subroutine finish_collisions

      use dissipation_coll_dougherty, only: finish_collisions_dougherty
      use dissipation_coll_fokkerplanck, only: finish_collisions_fp

      implicit none

      if (collision_model == "dougherty") then
         call finish_collisions_dougherty
      else if (collision_model == "fokker-planck") then
         call finish_collisions_fp
      end if

      initialised_collisions = .false.

   end subroutine finish_collisions

   !****************************************************************************
   subroutine advance_collisions_explicit(g, phi, bpar, gke_rhs)

      use mp, only: mp_abort
      use parameters_physics, only: full_flux_surface
      use stella_layouts, only: vmu_lo
      use grids_z, only: nzgrid
      use dissipation_coll_dougherty, only: advance_collisions_dougherty_explicit
      use dissipation_coll_fokkerplanck, only: advance_collisions_fp_explicit

      implicit none

      complex, dimension(:, :, -nzgrid:, :, vmu_lo%llim_proc:), intent(in) :: g
      complex, dimension(:, :, -nzgrid:, :), intent(in) :: phi, bpar
      complex, dimension(:, :, -nzgrid:, :, vmu_lo%llim_proc:), intent(in out) :: gke_rhs

      if (full_flux_surface) then
         call mp_abort("collisions not currently supported for full_flux_surface=T.  Aborting.")
      end if

      if (collision_model == "dougherty") then
         call advance_collisions_dougherty_explicit(g, phi, bpar, gke_rhs, time_collisions)
      else if (collision_model == "fokker-planck") then
         call advance_collisions_fp_explicit(g, phi, bpar, gke_rhs, time_collisions)
      end if

   end subroutine advance_collisions_explicit

   !****************************************************************************
   subroutine advance_collisions_implicit(mirror_implicit, phi, apar, bpar, g)

      use mp, only: proc0
      use redistribute, only: gather, scatter
      use calculations_redistribute, only: kxkyz2vmu
      use job_manage, only: time_message
      use grids_z, only: nzgrid
      use grids_velocity, only: set_vpa_weights
      use stella_layouts, only: vmu_lo
      use arrays_store_distribution_fn, only: gvmu
      use dissipation_coll_dougherty, only: advance_collisions_dougherty_implicit
      use dissipation_coll_fokkerplanck, only: advance_collisions_fp_implicit

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
