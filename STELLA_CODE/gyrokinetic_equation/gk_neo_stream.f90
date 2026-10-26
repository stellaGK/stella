! ================================================================================================================================================================================= !
! ----------------------------------------------------------- Evolves neoclassical corrections to parallel streaming. ------------------------------------------------------------- !​
! ================================================================================================================================================================================= !
! 
! This module evolves the following higher order neoclassical corrections: 
!
! =             
! 
! ================================================================================================================================================================================= !

module gk_neo_stream
   
   implicit none

   ! Make routines available to other modules. 
   public :: initialised_neo_stream
   public :: init_neo_stream
   public :: finish_neo_stream
   public :: advance_neo_stream_explicit

   private

   ! Only initialise once.
   logical :: initialised_neo_stream = .false.
   integer, dimension(:), allocatable :: neo_stream_sign
   integer, dimension(:), allocatable :: neo_stream_apar_sign

contains

! ================================================================================================================================================================================= !
! ------------------------------------------------------------- Initialise the neoclassical streaming correction. ----------------------------------------------------------------- ! 
! ================================================================================================================================================================================= !

    subroutine init_neo_stream
        ! Parallelisation.
        use mp, only: mp_abort
        use parallelisation_layouts, only: vmu_lo, iv_idx, imu_idx, is_idx

        ! Grids.
        use grids_time, only: code_dt
        use grids_species, only: spec, nspec
        use grids_velocity, only: maxwell_vpa, maxwell_mu, maxwell_fac
        use grids_velocity, only: vpa, mu, nvpa, nmu, vperp2
        use grids_z, only: nzgrid, nztot
        use grids_kxky, only: nalpha

        ! Geometry.
        use geometry, only: bmag, dbdzed, b_dot_gradz

        ! Parameters
        use parameters_physics, only: include_apar

        ! Arrays.
        use arrays, only: neo_stream, neo_stream_apar, initialised_neo_stream

        ! NEO data.
        use neoclassical_terms_neo, only: neo_vpa_fac_global

        ! For switching neoclassical streaming on and off.
        use parameters_physics, only: neostreamknob

        implicit none

        ! Local variables. 
        integer :: ia, iz, imu, iv, is

        ! Only intialise once.
        if (initialised_neo_stream) return
        initialised_neo_stream = .true.

        ! Allocate neo_stream = neo_stream[ialpha, iz, ivpa, imu, is].
        if (.not. allocated(neo_stream)) then
            allocate (neo_stream(nalpha, -nzgrid:nzgrid, nvpa, nmu, nspec)); neo_stream = 0.0
        end if

        ! Allocate neo_stream_apar = neo_stream_apar[ialpha, iz, ivpa, imu, is].
        if (.not. allocated(neo_stream_apar)) then
            allocate (neo_stream_apar(nalpha, -nzgrid:nzgrid, nvpa, nmu, nspec)); neo_stream_apar = 0.0
        end if

        ! Allocate neo_stream_sign = neo_stream_sign[i[vpa]]
        if (.not. allocated(neo_stream_sign)) then
            allocate (neo_stream_sign(nvpa)); neo_stream_sign = 0.0
        end if

        ! Allocate neo_stream_apar_sign = neo_stream_apar_sign[i[vpa]]
        if (.not. allocated(neo_stream_apar_sign)) then
            allocate (neo_stream_apar_sign(nvpa)); neo_stream_apar_sign = 0.0
        end if

        ! Calculate the higher order streaming coefficient for phi and bpar.
        do ia = 1, nalpha
            do iz = -nzgrid, nzgrid
                do iv = 1, nvpa
                    do imu = 1, nmu
                        do is = 1, nspec
                            neo_stream(ia, iz, iv, imu, is) = neostreamknob * code_dt &
                            * b_dot_gradz(ia, iz) * spec(is)%stm_psi0 * spec(is)%zt &
                            * 0.5 * neo_vpa_fac_global(iz, iv, imu, is, 1) &
                            * maxwell_fac(is) * maxwell_vpa(iv, is) &
                            * maxwell_mu(ia, iz, imu, is)
                        end do
                    end do
                end do
            end do
        end do

        if (include_apar) then
            ! Calculate the higher order streaming coefficient for apar.
            do ia = 1, nalpha
                do iz = -nzgrid, nzgrid
                    do iv = 1, nvpa
                        do imu = 1, nmu
                            do is = 1, nspec
                                neo_stream_apar(ia, iz, iv, imu, is) = - neostreamknob &
                                * code_dt &
                                * vpa(iv) * b_dot_gradz(ia, iz) &
                                * spec(is)%stm_psi0 * spec(is)%stm_psi0 * spec(is)%zt &
                                * neo_vpa_fac_global(iz, iv, imu, is, 1) &
                                * maxwell_fac(is) * maxwell_vpa(iv, is) &
                                * maxwell_mu(ia, iz, imu, is)
                            end do
                        end do
                    end do
                end do
            end do
        end if

        ! Calculate the sign of the streaming term for phi and bpar
        ! Do this for each point in the velocity space.
        do iv = 1, nvpa
            neo_stream_sign(iv) = int(sign(1.0, neo_stream(1, 0, iv, 1, 1))) 
        end do 

        if (include_apar) then
            ! Calculate the sign of the streaming term for apar.
            ! Do this for each point in the velocity space.
            do iv = 1, nvpa
                neo_stream_apar_sign(iv) = int(sign(1.0, neo_stream_apar(1, 0, iv, 1, 1))) 
            end do
        end if 

    end subroutine init_neo_stream

! ================================================================================================================================================================================= !
! ------------------------------------------------------------------------- Advance the terms explicitly. ------------------------------------------------------------------------- ! 
! ================================================================================================================================================================================= !

    subroutine advance_neo_stream_explicit(phi, apar, bpar, gout)
        ! Parallelisation.
        use mp, only: proc0, iproc, mp_abort, barrier
        use parallelisation_layouts, only: vmu_lo, iv_idx, imu_idx, is_idx

        ! Data arrays.
        use arrays, only: neo_stream, neo_stream_apar

        ! Grids.
        use grids_species, only: spec
        use grids_z, only: nzgrid, ntubes
        use grids_kxky, only: naky, nakx
        use grids_velocity, only: mu, vpa
        use grids_velocity, only: maxwell_vpa, maxwell_mu, maxwell_fac

        ! For calculating ∂<Χ_k>/∂z.
        use gk_parallel_streaming, only: get_dgdz_centered

        ! Parameters
        use parameters_physics, only: fphi, include_apar, include_bpar

        ! Calculations.
        use calculations_gyro_averages, only: gyro_average, gyro_average_j1
        use calculations_add_explicit_terms, only: add_explicit_term

        ! Time this routine.
        use timers, only: time_gke
        use job_manage, only: time_message

        implicit none

        complex, dimension(:, :, -nzgrid:, :), intent(in) :: phi, apar, bpar
        complex, dimension(:, :, -nzgrid:, :, vmu_lo%llim_proc:), intent(in out) :: gout

        ! Local variables.
        integer :: ivmu, iv, imu, is, ia, iz, it
        complex, dimension(:, :, :, :), allocatable :: g0, dgphi_dz, dgapar_dz, dgbpar_dz

        ! Allocate arrays needed for intermmediate calculations
        allocate (g0(naky, nakx, -nzgrid:nzgrid, ntubes))
        allocate (dgphi_dz(naky, nakx, -nzgrid:nzgrid, ntubes))
        allocate (dgapar_dz(naky, nakx, -nzgrid:nzgrid, ntubes))
        allocate (dgbpar_dz(naky, nakx, -nzgrid:nzgrid, ntubes))
        dgphi_dz = 0.0
        dgapar_dz = 0.0
        dgbpar_dz = 0.0

        ! Assume only one field line.
        ia = 1

        ! Start timing the advance.
        if (proc0) call time_message(.false., time_gke(:, 6), 'neo_stream advance')

        ! Calculate the parallel derivative of each field and add to the RHS of the GKE.
        do ivmu = vmu_lo%llim_proc, vmu_lo%ulim_proc
            iv = iv_idx(vmu_lo, ivmu)
            imu = imu_idx(vmu_lo, ivmu)
            is = is_idx(vmu_lo, ivmu)

            ! Gyroaverage.
            call gyro_average(phi, ivmu, g0(:, :, :, :))

            ! Get zed derivative.
            call get_dgdz_centered_neo(g0, ivmu, dgphi_dz)

            ! If apar is present, we calculate and add the correciton.
            if (include_apar) then
                ! Gyroaverage.
                call gyro_average(apar, ivmu, g0(:, :, :, :))

                ! Get the zed derivative.
                call get_dgdz_centered_neo_apar(g0, ivmu, dgapar_dz)
            end if

            ! If bpar is present, we must account for this too.
            if (include_bpar) then
                ! Gyroaverage.
                call gyro_average_j1(bpar, ivmu, g0(:, :, :, :))

                ! Get the zed derivative.
                call get_dgdz_centered_neo(g0, ivmu, dgbpar_dz)
            end if

            ! Construct the full parallel derivative.
            g0(:, :, :, :) = dgphi_dz(:, :, :, :) &
            + 4.0 * mu(imu) * spec(is)%tz * dgbpar_dz(:, :, :, :)

            ! Add neo_stream term to the right-hand-side of the GKE.
            call add_neo_stream_term(g0, ivmu, gout(:, :, :, :, ivmu))

            ! If we are including apar, add this term seperately.
            if (include_apar) then
                g0(:, :, :, :) = dgapar_dz(:, :, :, :)

                ! Add neo_stream term to the right-hand-side of the GKE.
                call add_neo_stream_apar_term(g0, ivmu, gout(:, :, :, :, ivmu))
            end if
        end do

        ! Deallocate temporary arrays.
        deallocate(g0)
        deallocate(dgphi_dz)
        deallocate(dgapar_dz)
        deallocate(dgbpar_dz)

        ! Stop timing the advance.
        if (proc0) call time_message(.false., time_gke(:, 6), 'neo_stream advance')
    end subroutine advance_neo_stream_explicit


! ================================================================================================================================================================================= !
! ------------------------------------------------------------------------------- Finish the terms. ------------------------------------------------------------------------------- ! 
! ================================================================================================================================================================================= !

    subroutine finish_neo_stream
        use arrays, only: neo_stream, neo_stream_apar, initialised_neo_stream

        implicit none

        if (allocated(neo_stream)) deallocate(neo_stream)
        if (allocated(neo_stream_apar)) deallocate(neo_stream_apar)
        if (allocated(neo_stream_sign)) deallocate(neo_stream_sign)
        initialised_neo_stream = .false.

    end subroutine finish_neo_stream


! ================================================================================================================================================================================= !
! --------------------------------------------------------------------------------------------------------------------------------------------------------------------------------- ! 
! ================================================================================================================================================================================= !

   ! Get second order accurate centered dg/dz, assuming delta zed is equally spaced
   subroutine get_dgdz_centered_neo(g, ivmu, dgdz)

      use calculations_finite_differences, only: second_order_centered_zed
      use parallelisation_layouts, only: vmu_lo
      use parallelisation_layouts, only: iv_idx
      use grids_z, only: nzgrid, delzed, ntubes
      use grids_extended_zgrid, only: neigen, nsegments
      use grids_extended_zgrid, only: iz_low, iz_up
      use grids_extended_zgrid, only: ikxmod
      use grids_extended_zgrid, only: fill_zed_ghost_zones
      use grids_extended_zgrid, only: periodic
      use grids_kxky, only: naky

      implicit none

      complex, dimension(:, :, -nzgrid:, :), intent(in) :: g
      complex, dimension(:, :, -nzgrid:, :), intent(out) :: dgdz
      integer, intent(in) :: ivmu

      integer :: iseg, ie, iky, iv, it
      complex, dimension(2) :: gleft, gright

      !-------------------------------------------------------------------------

      iv = iv_idx(vmu_lo, ivmu)
      do iky = 1, naky
         do it = 1, ntubes
            do ie = 1, neigen(iky)
               do iseg = 1, nsegments(ie, iky)
                  ! First fill in ghost zones at boundaries in g(z)
                  call fill_zed_ghost_zones(it, iseg, ie, iky, g(:, :, :, :), gleft, gright)
                  ! Now get dg/dz
                  call second_order_centered_zed(iz_low(iseg), iseg, nsegments(ie, iky), &
                     g(iky, ikxmod(iseg, ie, iky), iz_low(iseg):iz_up(iseg), it), &
                     delzed(0), neo_stream_sign(iv), gleft, gright, periodic(iky), &
                     dgdz(iky, ikxmod(iseg, ie, iky), iz_low(iseg):iz_up(iseg), it))
               end do
            end do
         end do
      end do

   end subroutine get_dgdz_centered_neo

! ================================================================================================================================================================================= !
! --------------------------------------------------------------------------------------------------------------------------------------------------------------------------------- ! 
! ================================================================================================================================================================================= !

   ! Get second order accurate centered dg/dz, assuming delta zed is equally spaced
   subroutine get_dgdz_centered_neo_apar(g, ivmu, dgdz)

      use calculations_finite_differences, only: second_order_centered_zed
      use parallelisation_layouts, only: vmu_lo
      use parallelisation_layouts, only: iv_idx
      use grids_z, only: nzgrid, delzed, ntubes
      use grids_extended_zgrid, only: neigen, nsegments
      use grids_extended_zgrid, only: iz_low, iz_up
      use grids_extended_zgrid, only: ikxmod
      use grids_extended_zgrid, only: fill_zed_ghost_zones
      use grids_extended_zgrid, only: periodic
      use grids_kxky, only: naky

      implicit none

      complex, dimension(:, :, -nzgrid:, :), intent(in) :: g
      complex, dimension(:, :, -nzgrid:, :), intent(out) :: dgdz
      integer, intent(in) :: ivmu

      integer :: iseg, ie, iky, iv, it
      complex, dimension(2) :: gleft, gright

      !-------------------------------------------------------------------------

      iv = iv_idx(vmu_lo, ivmu)
      do iky = 1, naky
         do it = 1, ntubes
            do ie = 1, neigen(iky)
               do iseg = 1, nsegments(ie, iky)
                  ! First fill in ghost zones at boundaries in g(z)
                  call fill_zed_ghost_zones(it, iseg, ie, iky, g(:, :, :, :), gleft, gright)
                  ! Now get dg/dz
                  call second_order_centered_zed(iz_low(iseg), iseg, nsegments(ie, iky), &
                     g(iky, ikxmod(iseg, ie, iky), iz_low(iseg):iz_up(iseg), it), &
                     delzed(0), neo_stream_apar_sign(iv), gleft, gright, periodic(iky), &
                     dgdz(iky, ikxmod(iseg, ie, iky), iz_low(iseg):iz_up(iseg), it))
               end do
            end do
         end do
      end do

   end subroutine get_dgdz_centered_neo_apar

! ==============================================================================================================================================$
! ----------------------------------------------------------------------------------------------------------------------------------------------$
! ==============================================================================================================================================$

   subroutine add_neo_stream_term(g, ivmu, src)
      ! Parallelisation.
      use parallelisation_layouts, only: vmu_lo
      use parallelisation_layouts, only: iv_idx, imu_idx, is_idx

      ! MP.
      use mp, only: mp_abort

      ! Grids.
      use grids_z, only: nzgrid

      ! Arrays.
      use arrays, only: neo_stream

      implicit none

      complex, dimension(:, :, -nzgrid:, :), intent(in) :: g
      complex, dimension(:, :, -nzgrid:, :), intent(in out) :: src
      integer, intent(in) :: ivmu

      integer :: iz, iv, is, imu

      ! ================================================================================ ! 

      iv = iv_idx(vmu_lo, ivmu)
      imu = imu_idx(vmu_lo, ivmu)
      is = is_idx(vmu_lo, ivmu)
      do iz = -nzgrid, nzgrid
         src(:, :, iz, :) = src(:, :, iz, :) + neo_stream(1, iz, iv, imu, is) * g(:, :, iz, :)
      end do

   end subroutine add_neo_stream_term

! ==============================================================================================================================================$
! ----------------------------------------------------------------------------------------------------------------------------------------------$
! ==============================================================================================================================================$

   subroutine add_neo_stream_apar_term(g, ivmu, src)
      ! Parallelisation.
      use parallelisation_layouts, only: vmu_lo
      use parallelisation_layouts, only: iv_idx, imu_idx, is_idx

      ! MP.
      use mp, only: mp_abort

      ! Grids.
      use grids_z, only: nzgrid

      ! Arrays.
      use arrays, only: neo_stream_apar

      implicit none

      complex, dimension(:, :, -nzgrid:, :), intent(in) :: g
      complex, dimension(:, :, -nzgrid:, :), intent(in out) :: src
      integer, intent(in) :: ivmu

      integer :: iz, iv, is, imu

      ! ================================================================================ ! 

      iv = iv_idx(vmu_lo, ivmu)
      imu = imu_idx(vmu_lo, ivmu)
      is = is_idx(vmu_lo, ivmu)
      do iz = -nzgrid, nzgrid
         src(:, :, iz, :) = src(:, :, iz, :) + neo_stream_apar(1, iz, iv, imu, is) * g(:, :, iz, :)
      end do
   end subroutine add_neo_stream_apar_term


! ================================================================================================================================================================================= !
! --------------------------------------------------------------------------------------------------------------------------------------------------------------------------------- ! 
! ================================================================================================================================================================================= !

end module gk_neo_stream

