! ================================================================================================================================================================================= !
! --------------------------------- Evolves neoclassical corrections proportional to the z derivative of the gyroaveraged generalised, ∂<Χ_k>/∂z. --------------------------------- !​
! ================================================================================================================================================================================= !
! 
! This module evolves the following higher order neoclassical corrections: 
!
! =             
!         
! Define the neoclassical ∂<Χ_k>/∂z coefficient as: 
! 
! <neo_dchidz_coeff> =
!
! This must be multiplied by ∂<Χ_k>/∂z and then added to the RHS of the GKE.
! 
! ================================================================================================================================================================================= !

module gk_neo_stream
   
   implicit none

   ! Make routines available to other modules. 
   public :: initialised_neo_stream
   public :: init_neo_stream, finish_neo_stream
   public :: advance_neo_stream_explicit

   private
   
   ! Only initialise once.
   logical :: initialised_neo_stream = .false.
   integer, dimension(:), allocatable :: neo_stream_sign

contains

! ================================================================================================================================================================================= !
! ------------------------------------------------------------------ Initialise the neoclassical ∂<Χ_k>/∂z terms. ----------------------------------------------------------------- ! 
! ================================================================================================================================================================================= !

    subroutine init_neo_stream
        ! Parallelisation.
        use mp, only: mp_abort
        use parallelisation_layouts, only: vmu_lo, iv_idx, imu_idx, is_idx

        ! Grids.
        use grids_time, only: code_dt
        use grids_species, only: spec, nspec
        use grids_velocity, only: maxwell_vpa, maxwell_mu, maxwell_fac
        use grids_velocity, only: mu, vperp2, vpa, nvpa
        use grids_z, only: nzgrid, nztot
        use grids_kxky, only: nalpha

        ! Geometry.
        use geometry, only: bmag, dbdzed, b_dot_gradz

        ! Arrays.
        use arrays, only: neo_stream, initialised_neo_stream

        ! NEO data.
        use neoclassical_terms_neo, only: neo_vpa_fac

        implicit none

        ! Local variables. 
        integer :: iz, iv, is, imu, ivmu

        ! Only intialise once.
        if (initialised_neo_stream) return
        initialised_neo_stream = .true.

        ! Allocate neo_stream = neo_stream[ialpha, iz, i[mu,vpa,s]].
        if (.not. allocated(neo_stream)) then
            allocate (neo_stream(nalpha, -nzgrid:nzgrid, vmu_lo%llim_proc:vmu_lo%ulim_alloc)); neo_stream = 0.0
        end if

        ! Allocate neo_stream_sign = neo_stream_sign[iv].
        if (.not. allocated(neo_stream_sign)) then
            allocate (neo_stream_sign(nvpa)); neo_stream_sign = 0.0
        end if

        ! Iterate over velocity space.
        do ivmu = vmu_lo%llim_proc, vmu_lo%ulim_proc
            is = is_idx(vmu_lo, ivmu)
            imu = imu_idx(vmu_lo, ivmu)
            iv = iv_idx(vmu_lo, ivmu)

            do iz = -nzgrid, nzgrid
                neo_stream(:, iz, ivmu) = 0.5 * code_dt * spec(is)%zt * spec(is)%stm * b_dot_gradz(:, iz) &
                * maxwell_vpa(iv, is) * maxwell_mu(:, iz, imu, is) * maxwell_fac(is) * neo_vpa_fac(iz, ivmu, 1)
            end do
        end do

        ! ================================================================================================================================================= !
        ! Create an array that gives the sign of the streaming coeffient.                                                                                   ! 
        ! This is needed becuase the upwinding factor is dependent on the direction of advection.                                                           !
        ! Here, neo_stream_sign set to +/- 1 depending on the sign of the parallel streaming term.                                                          !
        ! NB: stream_sign = -1 corresponds to positive advection velocity.                                                                                  !
        ! We only need to consider ia=1, iz=0 and is=1 because alpha, z and species dependencies do not lead to change in sign of the streaming pre-factor. !
        ! ================================================================================================================================================= !
 
        ! Iterate over velocity space.
        do ivmu = vmu_lo%llim_proc, vmu_lo%ulim_proc
            iv = iv_idx(vmu_lo, ivmu)

            neo_stream_sign(iv) = int(sign(1.0, neo_stream(1, 0, ivmu)))
        end do

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

        ! Parameters
        use parameters_physics, only: include_apar, include_bpar

        ! Calculations.
        use calculations_gyro_averages, only: gyro_average, gyro_average_j1
        use calculations_add_explicit_terms, only: add_explicit_term

        ! Time this routine.
        use timers, only: time_gke
        use job_manage, only: time_message

        implicit none

        complex, dimension(:, :, -nzgrid:, :), intent(in) :: phi, apar, bpar
        complex, dimension(:, :, -nzgrid:, :, vmu_lo%llim_proc:), intent(in out) :: gout
        complex, dimension(:, :, :, :), allocatable :: g0, dphi_dz, dapar_dz, dbpar_dz        
        complex, dimension(:, :, :, :, :), allocatable :: dg0_dz

        integer :: iv, is, imu, ivmu, ia, iz

        ! Allocate temporary arrays.
        allocate (g0(naky, nakx, -nzgrid:nzgrid, ntubes))
        allocate (dphi_dz(naky, nakx, -nzgrid:nzgrid, ntubes))
        allocate (dapar_dz(naky, nakx, -nzgrid:nzgrid, ntubes))
        allocate (dbpar_dz(naky, nakx, -nzgrid:nzgrid, ntubes))
        allocate (dg0_dz(naky, nakx, -nzgrid:nzgrid, ntubes, vmu_lo%llim_proc:vmu_lo%ulim_alloc))

        ! Assume we only have one field line.
        ia = 1

        ! ======================================================================================= ! 
        ! --------------------------------------------------------------------------------------- !
        ! ======================================================================================= ! 
        !                                                                                         !
        ! Calculate the parallel derivative of each field:                                        !
        !                                                                                         !
        ! <dphi_dz> = ∂<phi_k>/∂z, <dbpar_dz> = ∂<bpar_k>/∂z                                      !
        !                                                                                         !
        ! Then construct the parallel derivative of the total gyrokinetic potential.              !
        ! Mutlipy this by neo_stream and add to the right-hand-side of the GKE:                   !
        !                                                                                         ! 
        ! add_explicit_term(g0, neo_stream(1, :, :), gout)                                        !
        !                                                                                         ! 
        ! ======================================================================================= !
        ! --------------------------------------------------------------------------------------- ! 
        ! ======================================================================================= !

        ! Start timing the time advance.
        if (proc0) call time_message(.false., time_gke(:, 6), 'neo_stream advance')

        ! Calculate the parallel derivative depending on which fields are active. 
        do ivmu = vmu_lo%llim_proc, vmu_lo%ulim_proc
            iv = iv_idx(vmu_lo, ivmu)
            imu = imu_idx(vmu_lo, ivmu)
            is = is_idx(vmu_lo, ivmu)

            ! Get <phi>.         
            call gyro_average(phi, ivmu, g0(:, :, :, :))

            ! Get d<phi>/dz.
            call get_dgdz_centered_neo(g0, ivmu, dphi_dz)

            if (include_apar) then 
                ! Get <apar>.         
                call gyro_average(apar, ivmu, g0(:, :, :, :))

                ! Get d<apar>/dz.
                call get_dgdz_centered_neo(g0, ivmu, dapar_dz)
            else
                dapar_dz = 0.0
            end if

            if (include_bpar) then
                call gyro_average_j1(bpar, ivmu, g0(:, :, :, :))

                call get_dgdz_centered_neo(g0, ivmu, dbpar_dz)
            else
                dbpar_dz = 0.0
            end if 

            ! Construct the full z derivative. 
            dg0_dz(:, :, :, :, ivmu) = dg0_dz(:, :, :, :, ivmu) &
            + dphi_dz(:, :, :, :) - 2.0 * spec(is)%stm * vpa(iv) * dapar_dz(:, :, :, :)  + 4.0 * mu(imu) * spec(is)%tz * dbpar_dz(:, :, :, :)
        end do

        ! Add the term to the right-hand-side of the GKE. 
        call add_explicit_term(dg0_dz, neo_stream(1, :, :), gout)

        ! Allocate temporary arrays. 
        deallocate (g0)                          
        deallocate (dg0_dz)
        deallocate (dphi_dz)
        deallocate (dapar_dz)
        deallocate (dbpar_dz)

        ! Stop timing the time advance.
        if (proc0) call time_message(.false., time_gke(:, 6), 'neo_stream advance')

    end subroutine advance_neo_stream_explicit


! ================================================================================================================================================================================= !
! -------------------------- Get centered dg/dz using the neo_stream_sign: Get second order accurate centered dg/dz, assuming delta zed is equally spaced. ------------------------ ! 
! ================================================================================================================================================================================= !

    subroutine get_dgdz_centered_neo(g, ivmu, dgdz)
      ! Calculations. 
      use calculations_finite_differences, only: second_order_centered_zed

      ! Parallelisation. 
      use parallelisation_layouts, only: vmu_lo
      use parallelisation_layouts, only: iv_idx

      ! Grids. 
      use grids_z, only: nzgrid, delzed, ntubes
      use grids_extended_zgrid, only: neigen, nsegments
      use grids_extended_zgrid, only: iz_low, iz_up
      use grids_extended_zgrid, only: ikxmod
      use grids_extended_zgrid, only: fill_zed_ghost_zones
      use grids_extended_zgrid, only: periodic
      use grids_kxky, only: naky

      implicit none

      complex, dimension(:, :, -nzgrid:, :), intent(in) :: g
      complex, dimension(:, :, -nzgrid:, :), intent(in out) :: dgdz
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
! ------------------------------------------------------------------------------- Finish the terms. ------------------------------------------------------------------------------- ! 
! ================================================================================================================================================================================= !

    subroutine finish_neo_stream
        use arrays, only: neo_stream, initialised_neo_stream

        implicit none

        if (allocated(neo_stream)) deallocate (neo_stream)
        if (allocated(neo_stream_sign)) deallocate (neo_stream_sign)
        initialised_neo_stream = .false.

    end subroutine finish_neo_stream


! ================================================================================================================================================================================= !
! --------------------------------------------------------------------------------------------------------------------------------------------------------------------------------- ! 
! ================================================================================================================================================================================= !


end module gk_neo_stream
