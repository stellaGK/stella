 ! ================================================================================================================================================================================ !
! -------------------------------------------------------------- Evolves neoclassical corrections proportional to <A∥_k>. --------------------------------------------------------- !​
! ================================================================================================================================================================================= !
!                                                                                                                                                                                   !
! ================================================================================================================================================================================= !

module gk_neo_mirror
 
   implicit none

   ! Make routines available to other modules. 
   public :: init_neo_mirror, finish_neo_mirror
   public :: advance_neo_mirror_explicit

   public :: neo_mirror
   public :: neo_mirror_sign

   private

   integer, dimension(:, :), allocatable :: neo_mirror_sign
   real, dimension(:, :, :, :), allocatable :: neo_mirror


contains

! ================================================================================================================================================================================= !
! -------------------------------------------------------------------- Initialise the neoclassical mirror term. ------------------------------------------------------------------- ! 
! ================================================================================================================================================================================= !

    subroutine init_neo_mirror
        ! Parallelisation.
        use mp, only: mp_abort
        use parallelisation_layouts, only: vmu_lo, iv_idx, imu_idx, is_idx

        ! Grids.
        use grids_time, only: code_dt
        use grids_species, only: spec, nspec
        use grids_velocity, only: maxwell_vpa, maxwell_mu, maxwell_fac
        use grids_velocity, only: mu, vpa, nmu, nvpa, vperp2
        use grids_z, only: nzgrid, nztot
        use grids_kxky, only: nalpha

        ! Geometry.
        use geometry, only: bmag, dbdzed, b_dot_gradz

        ! Neoclassical.
        use neoclassical_terms_neo, only: neo_vpa_fac
        use neoclassical_terms_neo, only: dneo_phi_dz

        use arrays, only: neo_mirror_apar
        use arrays, only: initialised_neo_mirror

        ! For switching streaming on and off.
        use parameters_physics, only: neomirrorknob

        implicit none

        ! Local variables. 
        integer :: ia, iz, iv, is, imu, ivmu

        ! Only intialise once.
        if (initialised_neo_mirror) return
        initialised_neo_mirror = .true.

        ! Allocate arrays.
        if (.not. allocated(neo_mirror)) allocate (neo_mirror(nalpha, -nzgrid:nzgrid, nmu, nspec)); neo_mirror = 0.0
        if (.not. allocated(neo_mirror_sign)) allocate (neo_mirror_sign(nalpha, -nzgrid:nzgrid)); neo_mirror_sign = 0.0
        if (.not. allocated(neo_mirror_apar)) allocate (neo_mirror_apar(nalpha, -nzgrid:nzgrid, vmu_lo%llim_proc:vmu_lo%ulim_alloc)); neo_mirror_apar = 0.0

        ! We get a new mirror-like operator acting on the distribution.
        ! This is referred to as the "parallel-E" mirror term and is treated distinctly from the
        ! usual "grad-B" mirror term.  
        ! The neo_mirror array has a sign consistent with being on RHS of GKE.
        ! it is the factor multiplying dh/dvpa in the mirror term.
        ! For this term we choose to use h as the distribution rather than g.
        ! There appears to be no algebraic or numerical advantage to using g over h for this new term,
        ! using g seems to needlessly complicate things.  
        do imu = 1, nmu
            do ia = 1, nalpha
                do iz = -nzgrid, nzgrid
                    neo_mirror(ia, iz, imu, :) = 0.5 * neomirrorknob * code_dt * spec%z &
                    * spec%stm_psi0 * b_dot_gradz(ia, iz) * dneo_phi_dz(iz) 
                end do
            end do
        end do
       
        ! We need the sign of the neo_mirror array for calculating the dh/dvpa component of the parallel-E
        ! mirror term. 
        do ia = 1, nalpha
            do iz = -nzgrid, nzgrid
                neo_mirror_sign(ia, iz) = int(sign(1.0, neo_mirror(ia, iz, 1, 1)))
            end do
        end do

        ! The usual grad-B mirror term picks up a correction depending on the apar field only. 
        ! Iterate over velocity space.
        do ivmu = vmu_lo%llim_proc, vmu_lo%ulim_proc
            is = is_idx(vmu_lo, ivmu)
            imu = imu_idx(vmu_lo, ivmu)
            iv = iv_idx(vmu_lo, ivmu)

            do iz = -nzgrid, nzgrid
                neo_mirror_apar(:, iz, ivmu) = neomirrorknob * code_dt * spec(is)%z * mu(imu) &
                * b_dot_gradz(:, iz) * dbdzed(:, iz) * neo_vpa_fac(iz, ivmu, 1) * maxwell_vpa(iv, is) &
                * maxwell_mu(:, iz, imu, is) * maxwell_fac(is) / (spec(is)%mass * vpa(iv) )
            end do
        end do

    end subroutine init_neo_mirror

! ================================================================================================================================================================================= !
! ------------------------------------------------------------------------- Advance the terms explicitly. ------------------------------------------------------------------------- ! 
! ================================================================================================================================================================================= !

    subroutine advance_neo_mirror_explicit(g, phi, apar, bpar, gout)
        ! MP.
        use mp, only: mp_abort
        
        ! Parallelisation.
        use mp, only: proc0
        use parallelisation_layouts, only: vmu_lo
        use parallelisation_layouts, only: is_idx, iv_idx, imu_idx
        use parallelisation_layouts, only: fields_kxkyz
        use parallelisation_layouts, only: kxyz_lo, kxkyz_lo, vmu_lo

        ! Arrays.
        use arrays_distribution_function, only: gvmu
        use arrays, only: neo_mirror_apar

        ! Numerical parameters. 
        use parameters_numerical, only: mirror_implicit
        use parameters_numerical, only: neoclassical_mirror_implicit
        use parameters_physics, only: fphi
        use parameters_physics, only: include_apar

        ! Grids. 
        use grids_z, only: nzgrid, ntubes
        use grids_kxky, only: naky, nakx
        use grids_species, only: spec
        use grids_velocity, only: vpa, nvpa
        use grids_velocity, only: nmu
      
        ! Calculations.
        use calculations_add_explicit_terms, only: add_explicit_term
        use calculations_gyro_averages, only: gyro_average
        use calculations_tofrom_ghf, only: g_to_h

        ! Operations. 
        use redistribute, only: gather, scatter
        use initialise_redistribute, only: kxkyz2vmu, kxyz2vmu

        ! Time this routine.
        use timers, only: time_gke
        use job_manage, only: time_message

        implicit none

        complex, dimension(:, :, -nzgrid:, :, vmu_lo%llim_proc:), intent(in out) :: g
        complex, dimension(:, :, -nzgrid:, :), intent(in) :: phi, apar, bpar
        complex, dimension(:, :, -nzgrid:, :, vmu_lo%llim_proc:), intent(in out) :: gout        

        ! Local variables.
        integer :: ivmu, iv, is, imu
        complex, dimension(:, :, :), allocatable :: h0v
        complex, dimension(:, :, :, :, :), allocatable :: h0x
        complex, dimension(:, :, :, :), allocatable :: field
        complex, dimension(:, :, :, :, :), allocatable :: g0

        ! Start timing the advance.
        if (proc0) call time_message(.false., time_gke(:, 6), 'neo_mirror advance')

        ! Allocate temporary arrays.
        allocate(h0v(nvpa, nmu, kxkyz_lo%llim_proc:kxkyz_lo%ulim_alloc))
        allocate(h0x(naky, nakx, -nzgrid:nzgrid, ntubes, vmu_lo%llim_proc:vmu_lo%ulim_alloc))

        ! First we construct the parallel-E mirror term and add it to the RHS of the GKE.
        ! There is currently no implicit treatment for this term as there appears to be no need for one.
        ! However there is a logical flag anticipating that this may be implemented in the future.
        ! For now, we always use neoclassical_mirror_implicit = .False..
        if (.not. neoclassical_mirror_implicit) then  
            ! The operator acts on h, but the incoming distribution is g; transform g to h.
            call g_to_h(g, phi, apar, bpar, fphi)

            if (.not. fields_kxkyz) then
                call scatter(kxkyz2vmu, g, gvmu)
            end if
        
            ! Get dh/dvpa and store in h0v.
            h0v = gvmu
            call get_dhdvpa_explicit(h0v)

            ! Swap layouts so that (z,kx,ky) are local.
            call gather(kxkyz2vmu, h0v, h0x)
         
            ! Calculate the full parallel-E mirror term and add to the RHS of the GKE.
            call add_neo_parrallel_E_mirror_term(h0x, gout)
                  
            ! Now transform back to g.
            call g_to_h(g, phi, apar, bpar, -fphi)
        else
            call mp_abort('neoclassical_mirror_implicit = T is not implemented yet.')
        end if

        ! Deallocate temporary arrays.
        deallocate(h0v)
        deallocate(h0x)

        ! We now need to add the second term, proportional to the vpa derivative of apar. 
        ! We only need to do this if the grad-B mirror term is being evolved explicitly.
        if (.not. mirror_implicit .and. include_apar) then 
            ! Allocate temporary arrays.
            allocate(g0(naky, nakx, -nzgrid:nzgrid, ntubes, vmu_lo%llim_proc:vmu_lo%ulim_alloc))
            allocate(field(naky, nakx, -nzgrid:nzgrid, ntubes))
            
            ! Iterate over the (mu,vpa,s) points. 
            do ivmu = vmu_lo%llim_proc, vmu_lo%ulim_proc
                is = is_idx(vmu_lo, ivmu)
                iv = iv_idx(vmu_lo, ivmu)
                imu = imu_idx(vmu_lo, ivmu)

                ! Calculate the apar field. 
                ! All other factors are included in neo_mirror_apar. 
                field = apar

                ! Gyroaverage.
                call gyro_average(field, ivmu, g0(:, :, :, :, ivmu))
            end do

            ! Add this term to the right-hand-side of the GKE. 
            call add_explicit_term(g0, neo_mirror_apar(1, :, :), gout)

            ! Deallocate temporary arrays.
            deallocate(g0)
            deallocate(field)
        end if

        ! Stop timing the advance.
        if (proc0) call time_message(.false., time_gke(:, 6), 'neo_mirror advance')

    end subroutine advance_neo_mirror_explicit


! ================================================================================================================================================================================= !
! --------------------------------------------------------------------- Finish the mirror correction. ----------------------------------------------------------------------------- ! 
! ================================================================================================================================================================================= !

    subroutine finish_neo_mirror
        use arrays, only: neo_mirror_apar, initialised_neo_mirror

        implicit none

        if (allocated(neo_mirror_apar)) deallocate(neo_mirror_apar)
        initialised_neo_mirror = .false.

    end subroutine finish_neo_mirror


! ================================================================================================================================================================================= !
! --------------------------------------------------------------------------------------------------------------------------------------------------------------------------------- ! 
! ================================================================================================================================================================================= !

    subroutine get_dhdvpa_explicit(h)
      ! Parallelisation.
      use parallelisation_layouts, only: kxkyz_lo, iz_idx, is_idx
      
      ! Calculations.
      use calculations_finite_differences, only: third_order_upwind

      ! Grids. 
      use grids_velocity, only: nvpa, nmu, dvpa

      implicit none

      complex, dimension(:, :, kxkyz_lo%llim_proc:), intent(in out) :: h

      integer :: ikxkyz, imu, iz, is
      complex, dimension(:), allocatable :: tmp

      ! ======================================================== !

      ! Allocate temporary array.
      allocate (tmp(nvpa))

      ! Calculate the vpa derivative.
      ! Iterate over kxkyz.
      do ikxkyz = kxkyz_lo%llim_proc, kxkyz_lo%ulim_proc
         iz = iz_idx(kxkyz_lo, ikxkyz)
         is = is_idx(kxkyz_lo, ikxkyz)
         do imu = 1, nmu
            call third_order_upwind(1, h(:, imu, ikxkyz), dvpa, neo_mirror_sign(1, iz), tmp)
            h(:, imu, ikxkyz) = tmp
         end do
      end do

      ! Deallocate temporary arrays. 
      deallocate (tmp)

   end subroutine get_dhdvpa_explicit


! ================================================================================================================================================================================= !
! --------------------------------------------------------------------------------------------------------------------------------------------------------------------------------- ! 
! ================================================================================================================================================================================= !

    subroutine add_neo_parrallel_E_mirror_term(h, src)
      ! Parallelisation.
      use parallelisation_layouts, only: vmu_lo
      use parallelisation_layouts, only: imu_idx, is_idx

      ! Grids. 
      use grids_z, only: nzgrid, ntubes
      use grids_kxky, only: nakx

      implicit none

      complex, dimension(:, :, -nzgrid:, :, vmu_lo%llim_proc:), intent(in) :: h
      complex, dimension(:, :, -nzgrid:, :, vmu_lo%llim_proc:), intent(in out) :: src

      integer :: imu, is, ivmu
      integer :: it, iz, ikx

      ! ======================================================== !

      do ivmu = vmu_lo%llim_proc, vmu_lo%ulim_proc
         imu = imu_idx(vmu_lo, ivmu)
         is = is_idx(vmu_lo, ivmu)
         do it = 1, ntubes
            do iz = -nzgrid, nzgrid
               do ikx = 1, nakx
                  src(:, ikx, iz, it, ivmu) = src(:, ikx, iz, it, ivmu) &
                  + neo_mirror(1, iz, imu, is) * h(:, ikx, iz, it, ivmu)
               end do
            end do
         end do
      end do

   end subroutine add_neo_parrallel_E_mirror_term


! ================================================================================================================================================================================= !
! --------------------------------------------------------------------------------------------------------------------------------------------------------------------------------- ! 
! ================================================================================================================================================================================= !

end module gk_neo_mirror
