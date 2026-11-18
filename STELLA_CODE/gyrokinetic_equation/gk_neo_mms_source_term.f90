 ! ================================================================================================================================================================================= !
! -------------------------------------------------------------- Evolves a specified MMS source term. --------------------------------------------------------- !​
! ================================================================================================================================================================================= !
 
module gk_neo_mms_source_term
 
   implicit none

   ! Make routines available to other modules. 
   public :: initialised_neo_mms_source_term
   public :: init_neo_mms_source_term, finish_neo_mms_source_term
   public :: advance_neo_mms_source_term_explicit

   private
   
   ! Only initialise once.
   logical :: initialised_neo_mms_source_term = .false.

contains

! ================================================================================================================================================================================= !
! -------------------------------------------------------------------- Initialise the neoclassical MMS source term. --------------------------------------------------------------- ! 
! ================================================================================================================================================================================= !

    subroutine init_neo_mms_source_term
        ! Parallelisation.
        use mp, only: mp_abort
        use parallelisation_layouts, only: vmu_lo, iv_idx, imu_idx, is_idx

        ! Grids.
        use grids_time, only: code_dt, code_time
        use grids_species, only: spec, nspec
        use grids_velocity, only: maxwell_vpa, maxwell_mu, maxwell_fac
        use grids_velocity, only: mu, vpa, nmu, nvpa, vperp2
        use grids_z, only: zed, nzgrid
        use grids_kxky, only: akx, aky, nalpha, zed0, theta0

        ! Geometry.
        use geometry, only: bmag
        use geometry, only: B_times_kappa_dot_gradx, B_times_gradB_dot_gradx
        use geometry, only: b_dot_gradz

        ! Neoclassical.
        use neoclassical_terms_neo, only: neo_vpa_fac

        ! Arrays.
        use arrays, only: neo_mms_source_term
        use arrays, only: initialised_neo_mms_source_term

        ! Constants.
        use constants, only: zi

        implicit none

        ! Local variables. 
        integer :: iz, iv, is, imu, ivmu
        real    :: pi = acos(-1.0)

        ! Only intialise once.
        if (initialised_neo_mms_source_term) return
        initialised_neo_mms_source_term = .true.

        ! Allocate neo_mms_source_term = neo_mms_source_term[ialpha, iz, i[mu,vpa,s]].
        if (.not. allocated(neo_mms_source_term)) then
            allocate (neo_mms_source_term(nalpha, -nzgrid:nzgrid, vmu_lo%llim_proc:vmu_lo%ulim_alloc)); neo_mms_source_term = 0.0
        end if

        ! Iterate over velocity space.
        do ivmu = vmu_lo%llim_proc, vmu_lo%ulim_proc
            is = is_idx(vmu_lo, ivmu)
            imu = imu_idx(vmu_lo, ivmu)
            iv = iv_idx(vmu_lo, ivmu)

            ! This is currently setup for the neo_stream testing. 
            do iz = -nzgrid, nzgrid
                ! First calculate the magnetic geometry factor. 
                neo_mms_source_term(:, iz, ivmu) = 4.0 * spec(is)%stm &
                * b_dot_gradz(:, iz) * spec(is)%zt / 17.0 

                ! Multiply by the neoclassical factor.
                neo_mms_source_term(:, iz, ivmu) = neo_mms_source_term(:, iz, ivmu) &
                * ( vpa(iv) -  vpa(iv) *  vpa(iv) *  vpa(iv) ) 

                ! Multiply by maxwellian and other factors.
                neo_mms_source_term(:, iz, ivmu) = neo_mms_source_term(:, iz, ivmu) &
                 * maxwell_vpa(iv, is) * maxwell_mu(:, iz, imu, is) * maxwell_fac(is) 
                
                ! Multiply by the phi factor.
                neo_mms_source_term(:, iz, ivmu) = neo_mms_source_term(:, iz, ivmu) &
                * sin(2 * pi * code_time) * code_dt &
                * aky(1) * akx(1) * ( cos(0.5 * (zed(iz) - zed0(1,1))) ** 3 ) &
                * sin(0.5 * (zed(iz) - zed0(1,1))) 

                ! Add the g time derivative term term. 
                neo_mms_source_term(:, iz, ivmu) = neo_mms_source_term(:, iz, ivmu) &
                + 2.0 * pi * cos(2.0 * pi * code_time) * code_dt &
                * aky(1) * akx(1) * ( cos(0.5 * (zed(iz) - zed0(1,1))) ** 4 ) &
                * maxwell_vpa(iv, is) * maxwell_mu(:, iz, imu, is) * maxwell_fac(is)
            end do
        end do

    end subroutine init_neo_mms_source_term

! ================================================================================================================================================================================= !
! ------------------------------------------------------------------------- Advance the terms explicitly. ------------------------------------------------------------------------- ! 
! ================================================================================================================================================================================= !

    subroutine advance_neo_mms_source_term_explicit(gout)
        ! Parallelisation.
        use parallelisation_layouts, only: vmu_lo

        ! Data arrays.
        use arrays, only: neo_mms_source_term

        ! Grids. 
        use grids_z, only: nzgrid, ntubes
        use grids_kxky, only: naky, nakx

        ! Calculations.
        use calculations_add_explicit_terms, only: add_explicit_term

        implicit none

        complex, dimension(:, :, -nzgrid:, :, vmu_lo%llim_proc:), intent(in out) :: gout        

        ! Local variables.
        complex, dimension(:, :, :, :, :), allocatable :: g0

        ! Allocate temporary arrays.
        allocate(g0(naky, nakx, -nzgrid:nzgrid, ntubes, vmu_lo%llim_proc:vmu_lo%ulim_alloc))

        ! The prefactor on the MMS term is simply 1.  
        g0 = 1.0

        ! Add this term to the right-hand-side of the GKE.
        call add_explicit_term(g0, neo_mms_source_term(1, :, :), gout)

        ! Deallocate temporary arrays.
        deallocate(g0)

    end subroutine advance_neo_mms_source_term_explicit


! ================================================================================================================================================================================= !
! --------------------------------------------------------------------- Finish the mms source term. ------------------------------------------------------------------------------- ! 
! ================================================================================================================================================================================= !

    subroutine finish_neo_mms_source_term
        use arrays, only: neo_mms_source_term
        use arrays, only: initialised_neo_mms_source_term

        implicit none

        if (allocated(neo_mms_source_term)) deallocate (neo_mms_source_term)
        initialised_neo_mms_source_term = .false.

    end subroutine finish_neo_mms_source_term

! ================================================================================================================================================================================= !
! --------------------------------------------------------------------------------------------------------------------------------------------------------------------------------- ! 
! ================================================================================================================================================================================= !

end module gk_neo_mms_source_term
