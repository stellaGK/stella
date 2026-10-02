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

   integer, dimension(:), allocatable :: neo_stream_sign
   
   ! Only initialise once.
   logical :: initialised_neo_stream = .false.

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

        ! Arrays.
        use arrays, only: neo_stream, initialised_neo_stream

        ! NEO data.
        use neoclassical_terms_neo, only: neo_vpa_fac, neo_vpa_fac_global
        use neoclassical_terms_neo, only: distribute_vmus_over_procs

        ! For switching neoclassical streaming on and off.
        use parameters_physics, only: neostreamknob

        implicit none

        ! Local variables. 
        integer :: ia, iz, iv, is, imu, ivmu
        real, dimension(:, :, :, :, :), allocatable :: neo_stream_global

        ! Only intialise once.
        if (initialised_neo_stream) return
        initialised_neo_stream = .true.

        ! Allocate neo_stream = neo_stream_global[ialpha, iz, i[mu,vpa,s]].
        if (.not. allocated(neo_stream_global)) then
            allocate (neo_stream_global(nalpha, -nzgrid:nzgrid, nvpa, nmu, nspec)); neo_stream_global = 0.0
        end if

        ! Allocate neo_stream = neo_stream[ialpha, iz, i[mu,vpa,s]].
        if (.not. allocated(neo_stream)) then
            allocate (neo_stream(nalpha, -nzgrid:nzgrid, vmu_lo%llim_proc:vmu_lo%ulim_alloc)); neo_stream = 0.0
        end if
 
        ! Allocate neo_stream_sign = neo_stream_sign[i[vpa]]
        if (.not. allocated(neo_stream_sign)) then
            allocate (neo_stream_sign(nvpa)); neo_stream_sign = 0.0
        end if

        ! Calculate the higher order streaming coeffecient.
        do iz = -nzgrid, nzgrid
            do iv = 1, nvpa
                do imu = 1, nmu
                    do is = 1, nspec
                        neo_stream_global(:, iz, iv, imu, is) = neostreamknob * code_dt * 0.5 * spec(is)%zt * spec(is)%stm * b_dot_gradz(:, iz) &
                        * neo_vpa_fac_global(iz, iv, imu, is, 1) * maxwell_vpa(iv, is) * maxwell_mu(:, iz, imu, is) * maxwell_fac(is) 
                    end do
                end do
            end do
        end do

        ! Calculate the sign of the streaming term at each point in the velocity space.
        do iv = 1, nvpa
            neo_stream_sign(iv) = int(sign(1.0, neo_stream_global(1, 0, iv, 1, 1)))
        end do

        ! Distribute over velocity space. 
        do ia = 1, nalpha
            do iz = -nzgrid, nzgrid
                call distribute_vmus_over_procs(neo_stream_global(ia, iz, :, :, :), neo_stream(ia, iz, :))   
            end do
        end do

       ! Deallocate temporary arrays. 
       deallocate(neo_stream_global)

    end subroutine init_neo_stream

! ================================================================================================================================================================================= !
! ------------------------------------------------------------------------- Advance the terms explicitly. ------------------------------------------------------------------------- ! 
! ================================================================================================================================================================================= !

    subroutine advance_neo_stream_explicit(phi, apar, bpar, gout)
        ! Parallelisation.
        use mp, only: proc0
        use parallelisation_layouts, only: vmu_lo, iv_idx, imu_idx, is_idx
      
        ! Data arrays.
        use arrays, only: neo_stream

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
        integer :: iz
        integer :: iv, imu, is, ivmu
        complex, dimension(:, :, :, :), allocatable :: field
        complex, dimension(:, :, :, :), allocatable :: g0, dphi_dz, dapar_dz, dbpar_dz
        complex, dimension(:, :, :, :, :), allocatable :: dchi_dz

        ! Allocate temporary arrays.
        allocate(field(naky, nakx, -nzgrid:nzgrid, ntubes))
        allocate(g0(naky, nakx, -nzgrid:nzgrid, ntubes))
        allocate(dphi_dz(naky, nakx, -nzgrid:nzgrid, ntubes))
        allocate(dapar_dz(naky, nakx, -nzgrid:nzgrid, ntubes))
        allocate(dbpar_dz(naky, nakx, -nzgrid:nzgrid, ntubes))
        allocate(dchi_dz(naky, nakx, -nzgrid:nzgrid, ntubes, vmu_lo%llim_proc:vmu_lo%ulim_alloc))

        ! ======================================================================================= ! 
        ! --------------------------------------------------------------------------------------- !
        ! ======================================================================================= ! 
        !                                                                                         !
        ! Calculate the parallel derivative of each field:                                        !
        !                                                                                         !
        ! <dphi_dz> = ∂<phi_k>/∂z, <dapar_dz> = ∂<apar_k>/∂z, <dbpar_dz> = ∂<bpar_k>/∂z           !
        !                                                                                         !
        ! Then construct the parallel derivative of the total gyrokinetic potential.              !
        ! Mutlipy this by neo_stream and add to the right-hand-side of the GKE:                   !
        !                                                                                         ! 
        ! add_explicit_term(g0, neo_stream(1, :, :), gout)                                        !
        !                                                                                         ! 
        ! ======================================================================================= !
        ! --------------------------------------------------------------------------------------- ! 
        ! ======================================================================================= !

        ! Start timing the advance.
        if (proc0) call time_message(.false., time_gke(:, 6), 'neo_stream advance')

        ! Calculate the parallel derivative of each field and add to the RHS of the GKE. 
        do ivmu = vmu_lo%llim_proc, vmu_lo%ulim_proc
            iv = iv_idx(vmu_lo, ivmu)
            imu = imu_idx(vmu_lo, ivmu)
            is = is_idx(vmu_lo, ivmu)

            ! Calculate phi.
            field = fphi * phi

            ! Gyroaverage.
            call gyro_average(field, ivmu, g0(:, :, :, :))

            ! Get the z derivative.
            call get_dgdz_centered_neo(g0, ivmu, dphi_dz)
        end do

        ! If apar is present, we calculate and add the correciton.
        if (include_apar) then
            do ivmu = vmu_lo%llim_proc, vmu_lo%ulim_proc
                iv = iv_idx(vmu_lo, ivmu)
                imu = imu_idx(vmu_lo, ivmu)
                is = is_idx(vmu_lo, ivmu)

                field = 2.0 * vpa(iv) * spec(is)%stm_psi0 * apar

                ! Gyroaverage.
                call gyro_average(field, ivmu, g0(:, :, :, :))

                ! Get the z derivative.
                call get_dgdz_centered_neo(g0, ivmu, dapar_dz)
            end do
        end if

        ! If bpar is present, we must account for this too.
        if (include_bpar) then
            do ivmu = vmu_lo%llim_proc, vmu_lo%ulim_proc
                iv = iv_idx(vmu_lo, ivmu)
                imu = imu_idx(vmu_lo, ivmu)
                is = is_idx(vmu_lo, ivmu)

                field = 4.0 * mu(imu) * spec(is)%tz * bpar
               
                ! Gyroaverage.
                call gyro_average_j1(field, ivmu, g0(:, :, :, :))

                ! Get the z derivative.
                call get_dgdz_centered_neo(g0, ivmu, dbpar_dz)
            end do
        end if

        ! Construct dchidz. 
        do ivmu = vmu_lo%llim_proc, vmu_lo%ulim_proc
            dchi_dz(:, :, :, :, ivmu) = dphi_dz 

            if (include_apar) then
                dchi_dz(:, :, :, :, ivmu) = dchi_dz(:, :, :, :, ivmu) + dapar_dz
            end if

            if (include_bpar) then
                dchi_dz(:, :, :, :, ivmu) = dchi_dz(:, :, :, :, ivmu) + dbpar_dz
            end if
        end do

        ! Add the term to the right-hand-side of the GKE. 
        call add_explicit_term(dchi_dz, neo_stream(1, :, :), gout)

        ! Deallocate temporary arrays.
        deallocate(field)
        deallocate(g0)                          
        deallocate(dphi_dz)
        deallocate(dapar_dz)
        deallocate(dbpar_dz)
        deallocate(dchi_dz)

        ! Stop timing the advance.
        if (proc0) call time_message(.false., time_gke(:, 6), 'neo_stream advance')

    end subroutine advance_neo_stream_explicit


! ================================================================================================================================================================================= !
! ------------------------------------------------------------------------------- Finish the terms. ------------------------------------------------------------------------------- ! 
! ================================================================================================================================================================================= !

    subroutine finish_neo_stream
        use arrays, only: neo_stream, initialised_neo_stream

        implicit none

        if (allocated(neo_stream)) deallocate(neo_stream)
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
                     delzed(0), neo_stream_sign(ivmu), gleft, gright, periodic(iky), &
                     dgdz(iky, ikxmod(iseg, ie, iky), iz_low(iseg):iz_up(iseg), it))
               end do
            end do
         end do
      end do

   end subroutine get_dgdz_centered_neo

! ================================================================================================================================================================================= !
! --------------------------------------------------------------------------------------------------------------------------------------------------------------------------------- ! 
! ================================================================================================================================================================================= !


end module gk_neo_stream
