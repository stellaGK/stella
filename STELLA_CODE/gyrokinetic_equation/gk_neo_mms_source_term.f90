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
   !   S(t) = [ sin(2 pi t / tau) * neo_mms_source_term + cos(2 pi t / tau) * mms_source_dgdt ] * code_dt
   ! neo_mms_source_term holds the part that cancels the neo streaming correction,
   ! mms_source_dgdt the part from dg_MS/dt.
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
        real    :: f, dfdz

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

        ! Manufactured solution for the neo_stream test (F1 = vpa^2, J0 = 0.5, adiabatic electrons, n = T = Z_i = 1):
        !   g_MS  = sin(2 pi t / tau) kx ky f(z) F0,
        !   Phi'  = (4/17) sin(2 pi t / tau) kx ky f(z),
        !   f(z)  = cos^4(z/2) [1 + sin(z)].
        ! cos^4(z/2) makes f vanish (with its first three derivatives) at z = +/- pi; the [1 + sin(z)]
        ! factor breaks the z -> -z symmetry so the field-line-averaged |g|^2 sees the streaming error.
        ! No zed0 shift. The source S = dg_MS/dt - (neo streaming correction of g_MS) is split as
        !   S = sin(2 pi t / tau) * neo_mms_source_term + cos(2 pi t / tau) * mms_source_dgdt,
        ! with the time factors applied in advance_neo_mms_source_term_explicit.
        do ivmu = vmu_lo%llim_proc, vmu_lo%ulim_proc
            is = is_idx(vmu_lo, ivmu)
            imu = imu_idx(vmu_lo, ivmu)
            iv = iv_idx(vmu_lo, ivmu)

            do iz = -nzgrid, nzgrid
                ! Profile and its z-derivative.
                f    = cos(0.5 * zed(iz))**4 * (1.0 + sin(zed(iz)))
                dfdz = -2.0 * cos(0.5 * zed(iz))**3 * sin(0.5 * zed(iz)) * (1.0 + sin(zed(iz))) &
                       + cos(0.5 * zed(iz))**4 * cos(zed(iz))

                ! Part cancelling the neo streaming correction (multiplies sin(2 pi t / tau)):
                ! the correction of g_MS is (2/17) vth b.grad(z) (Z/T) (vpa - vpa^3) F0 kx ky f'(z),
                ! so the source carries minus that.
                neo_mms_source_term(:, iz, ivmu) = -2.0 * spec(is)%stm &
                    * b_dot_gradz(:, iz) * spec(is)%zt / 17.0 &
                    * ( vpa(iv) - vpa(iv) * vpa(iv) * vpa(iv) ) &
                    * maxwell_vpa(iv, is) * maxwell_mu(:, iz, imu, is) * maxwell_fac(is) &
                    * aky(1) * akx(1) * dfdz

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
        ! (add_explicit_gyrokinetic_terms returns rhs*dt).
        source_now = ( sin(2.0 * pi * mms_stage_time / mms_tau) * neo_mms_source_term(1, :, :) &
                     + cos(2.0 * pi * mms_stage_time / mms_tau) * mms_source_dgdt(1, :, :) ) * code_dt

        ! The prefactor on the MMS term is simply 1.  
        g0 = 1.0

        ! Add this term to the right-hand-side of the GKE.
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