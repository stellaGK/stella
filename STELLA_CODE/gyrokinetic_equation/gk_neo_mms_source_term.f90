! ================================================================================================================================================================================= !
! -------------------------------------------------------------- Evolves a specified MMS source term. --------------------------------------------------------- !​
! ================================================================================================================================================================================= !
 
module gk_neo_mms_source_term
 
   implicit none

   ! Make routines available to other modules. 
   public :: initialised_neo_mms_source_term
   public :: init_neo_mms_source_term, finish_neo_mms_source_term
   public :: advance_neo_mms_source_term_explicit
   public :: mms_stage_time

   private
   
   ! Only initialise once.
   logical :: initialised_neo_mms_source_term = .false.

   ! Time at which the source is evaluated. Set by the explicit RK drivers before each stage
   ! (code_time is only advanced once per full step, so it cannot be used inside the stages).
   real :: mms_stage_time = 0.0

   ! Period of the MMS time dependence: g_MS ~ sin(2 pi t / mms_tau). A long period keeps the RK3 time
   ! error of g_MS (~ (delt/tau)^2) well below the spatial error of the term being tested.
   real, parameter :: mms_tau = 20.0

   ! Time-independent (z, velocity) factors of the source. The time dependence is applied in
   ! advance_neo_mms_source_term_explicit, at the stage time, every time the RHS is built:
   !   S(t) = [ -i sin(2 pi t / tau) * neo_mms_source_term + cos(2 pi t / tau) * mms_source_dgdt ] * code_dt
   ! neo_mms_source_term holds the (real) part that cancels the neo wdrifty correction (the factor -i
   ! is applied in the advance routine), mms_source_dgdt the part from dg_MS/dt.
   real, dimension(:, :, :), allocatable :: mms_source_dgdt

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
        real    :: f

        ! Only intialise once.
        if (initialised_neo_mms_source_term) return
        initialised_neo_mms_source_term = .true.

        ! Allocate neo_mms_source_term = neo_mms_source_term[ialpha, iz, i[mu,vpa,s]].
        if (.not. allocated(neo_mms_source_term)) then
            allocate (neo_mms_source_term(nalpha, -nzgrid:nzgrid, vmu_lo%llim_proc:vmu_lo%ulim_alloc)); neo_mms_source_term = 0.0
        end if

        if (.not. allocated(mms_source_dgdt)) then
            allocate (mms_source_dgdt(nalpha, -nzgrid:nzgrid, vmu_lo%llim_proc:vmu_lo%ulim_alloc)); mms_source_dgdt = 0.0
        end if

        ! Manufactured solution for the neo_wdrifty test (F1 = vpa^2, J0 = 0.5, adiabatic electrons, n = T = Z_i = 1):
        !   g_MS  = sin(2 pi t / tau) kx ky f(z) F0,
        !   Phi'  = (4/17) sin(2 pi t / tau) kx ky f(z),
        !   f(z)  = cos^4(z/2) [1 + sin(z)].
        ! (A global phase e^{i pi/4} on g_MS and Phi' changes nothing for a single ky > 0 mode, so it is omitted.)
        ! The neo wdrifty correction of g_MS is
        !   i ky (1/B) [b x (vpa^2 kappa + mu grad B)].grad(y) J0 Phi' F0 (1 - vpa^2)
        ! = i sin(2 pi t / tau) * (2/17) ky (1 - vpa^2) [vpa^2 B_times_kappa_dot_grady + mu B B_times_gradB_dot_grady] F0 kx ky f(z),
        ! using (1 - vpa^2) = d_vpa F1 shorthand for F1 = vpa^2, and stella's normalisations
        ! B_times_kappa_dot_grady = (b x kappa).grad(y)/B and B_times_gradB_dot_grady = (b x grad B).grad(y)/B^2.
        ! The source S = dg_MS/dt - (correction) is split as
        !   S = -i sin(2 pi t / tau) * neo_mms_source_term + cos(2 pi t / tau) * mms_source_dgdt,
        ! with the time factors and the -i applied in advance_neo_mms_source_term_explicit.
        do ivmu = vmu_lo%llim_proc, vmu_lo%ulim_proc
            is = is_idx(vmu_lo, ivmu)
            imu = imu_idx(vmu_lo, ivmu)
            iv = iv_idx(vmu_lo, ivmu)

            do iz = -nzgrid, nzgrid
                ! Profile.
                f = cos(0.5 * zed(iz))**4 * (1.0 + sin(zed(iz)))

                ! Real factor of the part cancelling the neo wdrifty correction (multiplies -i sin(2 pi t / tau)).
                neo_mms_source_term(:, iz, ivmu) = 2.0 / 17.0 * akx(1) &
                    * ( 1.0 - vpa(iv) * vpa(iv) ) &
                    * ( vpa(iv) * vpa(iv) * B_times_kappa_dot_gradx(:, iz) &
                      + mu(imu) * bmag(:, iz) * B_times_gradB_dot_gradx(:, iz) ) &
                    * maxwell_vpa(iv, is) * maxwell_mu(:, iz, imu, is) * maxwell_fac(is) &
                    * aky(1) * akx(1) * f

                ! Part from dg_MS/dt (multiplies cos(2 pi t / tau)): (2 pi / tau) kx ky f(z) F0.
                mms_source_dgdt(:, iz, ivmu) = 2.0 * pi / mms_tau &
                    * aky(1) * akx(1) * f &
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

        ! Time.
        use grids_time, only: code_dt

        ! Constants.
        use constants, only: zi

        implicit none

        complex, dimension(:, :, -nzgrid:, :, vmu_lo%llim_proc:), intent(in out) :: gout        

        ! Local variables.
        complex, dimension(:, :, :, :, :), allocatable :: g0
        real, dimension(:, :), allocatable :: source_now
        real :: pi

        pi = acos(-1.0)

        ! Allocate temporary arrays.
        allocate(g0(naky, nakx, -nzgrid:nzgrid, ntubes, vmu_lo%llim_proc:vmu_lo%ulim_alloc))
        allocate(source_now(-nzgrid:nzgrid, vmu_lo%llim_proc:vmu_lo%ulim_alloc))

        ! Source at the current RK stage time, times the current time step
        ! (add_explicit_gyrokinetic_terms returns rhs*dt). The source is complex, so it is added in two
        ! parts, with the complex constant carried by g0 and the real (z, velocity) factor by source_now.

        ! dg_MS/dt part: real, prefactor 1.
        source_now = cos(2.0 * pi * mms_stage_time / mms_tau) * mms_source_dgdt(1, :, :) * code_dt
        g0 = 1.0
        call add_explicit_term(g0, source_now, gout)

        ! Part cancelling the neo wdrifty correction: prefactor -i.
        source_now = sin(2.0 * pi * mms_stage_time / mms_tau) * neo_mms_source_term(1, :, :) * code_dt
        g0 = -zi
        call add_explicit_term(g0, source_now, gout)

        ! Deallocate temporary arrays.
        deallocate(g0)
        deallocate(source_now)

    end subroutine advance_neo_mms_source_term_explicit


! ================================================================================================================================================================================= !
! --------------------------------------------------------------------- Finish the mms source term. ------------------------------------------------------------------------------- ! 
! ================================================================================================================================================================================= !

    subroutine finish_neo_mms_source_term
        use arrays, only: neo_mms_source_term
        use arrays, only: initialised_neo_mms_source_term

        implicit none

        if (allocated(neo_mms_source_term)) deallocate (neo_mms_source_term)
        if (allocated(mms_source_dgdt)) deallocate (mms_source_dgdt)
        initialised_neo_mms_source_term = .false.

    end subroutine finish_neo_mms_source_term

! ================================================================================================================================================================================= !
! --------------------------------------------------------------------------------------------------------------------------------------------------------------------------------- ! 
! ================================================================================================================================================================================= !

end module gk_neo_mms_source_term