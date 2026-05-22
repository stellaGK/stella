! ================================================================================================================================================================================== !
! ------------------------------------- The following routines are responsible for calculating and advancing the field equations in the HO theory ---------------------------------- ! 
! -------------------------------------------------------- See Chapter 7 of: https://www.overleaf.com/read/kwrvhpvypbxq#322021. ---------------------------------------------------- !
! ================================================================================================================================================================================== !

module field_equations_fluxtube_neoclassical   
    implicit none

    ! Make routines available to other modules.
    public :: init_neo_electrostatic_fields
    public :: init_neo_electromagnetic_fields
    public :: finish_neo_electrostatic_fields
    public :: finish_neo_electromagnetic_fields
    public :: advance_apar_neo   

    public :: advance_fields_fluxtube_using_neo_field_equations

   ! Advance fields for g(kx,ky,z,ivpamus) and g(vpa,mu,ikxkyzs), depending on the layout.
   interface advance_fields_fluxtube_using_neo_field_equations
      module procedure advance_fields_fluxtube_using_neo_field_equations_vmulo
      module procedure advance_fields_fluxtube_using_neo_field_equations_kxkyzlo
   end interface

   private

contains


! ================================================================================================================================================================================ !
! ------------------------------------------------------ Advances the field equations based on the g(kx,ky,z,ivpamus) layout. ---------------------------------------------------- !
! ================================================================================================================================================================================ !

   subroutine advance_fields_fluxtube_using_neo_field_equations_vmulo(g, phi, apar, bpar, dist, skip_fsa)
      ! Parallelisation.
      use mp, only: mp_abort, proc0
      use job_manage, only: time_message
      use parallelisation_layouts, only: vmu_lo, iv_idx, imu_idx
      use timers, only: time_field_solve
      
      ! Arrays.
      use arrays_distribution_function, only: phi_gyro
      
      ! Parameters.
      use parameters_physics, only: include_apar, include_bpar
      use parameters_physics, only: fphi
      use parameters_physics, only: beta
      
      ! Grids.
      use grids_z, only: nzgrid
      use grids_species, only: spec
      use grids_velocity, only: vpa, mu      

      ! Calculations.
      use calculations_gyro_averages, only: gyro_average, gyro_average_j1
      use calculations_velocity_integrals, only: integrate_species
      
      implicit none

      ! Arguments
      complex, dimension(:, :, -nzgrid:, :, vmu_lo%llim_proc:), intent(in) :: g
      complex, dimension(:, :, -nzgrid:, :), intent(out) :: phi, apar, bpar
      character(*), intent(in) :: dist 
      logical, optional, intent(in) :: skip_fsa
      
      ! Local variables.
      logical :: skip_fsa_local
      integer :: iv, ivmu, imu
      
      ! =============================================================================== !
      
      ! Used for the Dougherty collision operator.
      skip_fsa_local = .false.
      if (present(skip_fsa)) skip_fsa_local = skip_fsa

      ! Initialise the electrostatic potential phi.
      phi = 0.
      
      if (fphi > epsilon(0.0)) then
          ! Start timer.
          if (proc0) call time_message(.false., time_field_solve(:, 3), ' int_dv_g')
         
          ! First gyro-average the distribution function g at each phase space location and store this as phi_gyro = <g>_R = J_0 g in k-space.
          call gyro_average(g, phi_gyro)
        
          ! Calculate phi = sum_s Z_s n_s [ (2B/sqrt(pi)) int dvpa int dmu J_0 * g ] = integrate_species( J_0 * g ).
          call integrate_species(phi_gyro, spec%z * spec%dens_psi0, phi)

          ! Stop timer.
          if (proc0) call time_message(.false., time_field_solve(:, 3), ' int_dv_g')         
      end if

      ! Initialise the parallel magnetic potential apar.
      apar = 0.

      if (include_apar) then
          ! If fphi > 0, then phi_gyro = <g> already calculated above.
          call gyro_average(g, phi_gyro)
         
          ! For parallel Amperes Law, need to calculate parallel current rather than density, so multiply <g> by vpa before integrating. 
          ! Because we are parallelising over (vpa, mu) we need to first get the vpa index, then multiply with vpa.
          do ivmu = vmu_lo%llim_proc, vmu_lo%ulim_proc
              iv = iv_idx(vmu_lo, ivmu)
              phi_gyro(:, :, :, :, ivmu) = phi_gyro(:, :, :, :, ivmu) * vpa(iv)
          end do
         
          ! Integrate vpa*<g> over velocity space and sum over species store result in apar, which will be further modified below to account for field coupling.
          call integrate_species(phi_gyro, beta * spec%z * spec%dens_psi0 * spec%stm_psi0, apar)
         
          ! End timer.
          if (proc0) call time_message(.false., time_field_solve(:, 3), ' int_dv_g')       
      end if

      ! Initialise the parallel magnetic field bpar.
      bpar = 0.

      if (include_bpar) then
          ! Gyroaverage the distribution function g at each phase space location
          call gyro_average_j1(g, phi_gyro)

          ! Multiply by mu factor from vperp2.
          do ivmu = vmu_lo%llim_proc, vmu_lo%ulim_proc
              imu = imu_idx(vmu_lo, ivmu)
              phi_gyro(:, :, :, :, ivmu) = phi_gyro(:, :, :, :, ivmu) * mu(imu)
          end do
         
          ! Integrate <g> over velocity space and sum over species and store result in bpar, which will be further modified below to account for field coupling.
          call integrate_species(phi_gyro, -2.0 * beta * spec%temp_psi0 * spec%dens_psi0, bpar)
      end if

      ! ========================================================================================================================================================= !
      ! ---------- Now that we have the moments of gneo, we need to get the correct fields. This will depend on which fields we are running with. --------------- ! 
      ! ========================================================================================================================================================= !

      ! If we are running electrostatic, we only need to get phi and there is not coupling between fields. Calculate phi by dividing with denominator_fields_neo[iky,ikz,iz].  
      if (fphi > epsilon(0.0) .and. .not. include_apar .and. .not. include_bpar) then     
          call calculate_neo_phi(phi, dist, skip_fsa_local)
      end if  

      ! If we are running electromagnetic but no bpar, the phi and apar field equations are a 2 x 2 matrix solve.
      if (fphi > epsilon(0.0) .and. include_apar .and. .not. include_bpar) then
          call calculate_neo_phi_and_apar(phi, apar, dist, skip_fsa_local)
      end if

      ! If we are running electromagnetic fully electromagnetic, the phi apar and bpar field equations are a 3 x 3 matrix. 
      if (fphi > epsilon(0.0) .and. include_apar .and. include_bpar) then
          call calculate_neo_phi_apar_and_bpar(phi, apar, bpar, dist, skip_fsa_local)
      end if 

   end subroutine advance_fields_fluxtube_using_neo_field_equations_vmulo


! ================================================================================================================================================================================ !
! ------------------------------------------------------ Advances the field equations based on the g(kx,ky,z,ikxkyzs) layout. ---------------------------------------------------- !
! ================================================================================================================================================================================ !

   subroutine advance_fields_fluxtube_using_neo_field_equations_kxkyzlo(g, phi, apar, bpar, dist, skip_fsa)
      ! Parallelisation.
      use mp, only: proc0
      use mp, only: sum_allreduce, mp_abort
      use job_manage, only: time_message
      use parallelisation_layouts, only: kxkyz_lo
      use parallelisation_layouts, only: iz_idx, it_idx, ikx_idx, iky_idx, is_idx
      use timers, only: time_field_solve
      
      ! Parameters.
      use parameters_physics, only: fphi, beta
      use parameters_physics, only: include_apar, include_bpar
      use parameters_physics, only: radial_variation
      
      ! Grids.
      use grids_velocity, only: vpa, mu
      use grids_velocity, only: nvpa, nmu
      use grids_species, only: spec
      use grids_z, only: nzgrid
      
      ! Calculations.
      use calculations_gyro_averages, only: gyro_average, gyro_average_j1
      use calculations_velocity_integrals, only: integrate_vmu

      implicit none

      ! Arguments.
      complex, dimension(:, :, kxkyz_lo%llim_proc:), intent(in) :: g
      complex, dimension(:, :, -nzgrid:, :), intent(inout) :: phi, apar, bpar
      logical, optional, intent(in) :: skip_fsa
      character(*), intent(in) :: dist
      
      ! Local variables.
      complex :: tmp
      real :: wgt
      complex, dimension(:, :), allocatable :: g0
      integer :: ikxkyz, iz, it, ikx, iky, is, ia
      logical :: skip_fsa_local
      
      ! =============================================================================== !
      
      ! Used for the Dougherty collision operator.
      skip_fsa_local = .false.
      if (present(skip_fsa)) skip_fsa_local = skip_fsa

      ! Assume we only have one field line.
      ia = 1
      
      ! Initialise the electrostatic potential phi.
      phi = 0.

      if (fphi > epsilon(0.0)) then
          ! Start timer.
          if (proc0) call time_message(.false., time_field_solve(:, 3), ' int_dv_g')
         
          ! Allocate temporary arrays.
          allocate (g0(nvpa, nmu))
         
          ! Iterate over the (kx,ky,z,mu,vpa,s) points.
          do ikxkyz = kxkyz_lo%llim_proc, kxkyz_lo%ulim_proc
              iz = iz_idx(kxkyz_lo, ikxkyz)
              it = it_idx(kxkyz_lo, ikxkyz)
              ikx = ikx_idx(kxkyz_lo, ikxkyz)
              iky = iky_idx(kxkyz_lo, ikxkyz)
              is = is_idx(kxkyz_lo, ikxkyz)
            
              ! First gyro-average the distribution function g at each phase space location and store this as g0 = <g>_R = J_0 g.
              call gyro_average(g(:, :, ikxkyz), ikxkyz, g0)
            
              ! Calculate phi = sum_s Z_s n_s [ (2B/sqrt(pi)) int dvpa int dmu J_0 * g ].
              wgt = spec(is)%z * spec(is)%dens_psi0
              call integrate_vmu(g0, iz, tmp)
              phi(iky, ikx, iz, it) = phi(iky, ikx, iz, it) + wgt * tmp           
          end do
         
          ! Deallocate temporary array.
          deallocate (g0)
         
          ! Sum the values on all processors and send them to <proc0>.
          call sum_allreduce(phi)         
      end if

      ! Initialise the apar field
      apar = 0.
      
      ! Calculate the apar moment. 
      if (include_apar) then
          ! Allocate temporary array.
          allocate (g0(nvpa, nmu))
         
          ! Iterate over the (kx,ky,z,mu,vpa,s) points
          ! This gives: 2 β sum_s Z_s n_s vth J_0 g 
          do ikxkyz = kxkyz_lo%llim_proc, kxkyz_lo%ulim_proc
              iz = iz_idx(kxkyz_lo, ikxkyz)
              it = it_idx(kxkyz_lo, ikxkyz)
              ikx = ikx_idx(kxkyz_lo, ikxkyz)
              iky = iky_idx(kxkyz_lo, ikxkyz)
              is = is_idx(kxkyz_lo, ikxkyz)
            
              call gyro_average(spread(vpa, 2, nmu) * g(:, :, ikxkyz), ikxkyz, g0)
              wgt = 2.0 * beta * spec(is)%z * spec(is)%dens * spec(is)%stm
              call integrate_vmu(g0, iz, tmp)
              apar(iky, ikx, iz, it) = apar(iky, ikx, iz, it) + tmp * wgt
          end do
         
          ! Sum the values on all processors and send them to <proc0>.
          call sum_allreduce(apar)
          deallocate (g0)
      end if

      ! Initialise the bpar field
      bpar = 0.

      ! Calculate the bpar moment. 
      if (include_bpar) then
          ! Allocate temporary array.
          allocate (g0(nvpa, nmu))
          
          ! Iterate over the (kx,ky,z,mu,vpa,s) points.
          do ikxkyz = kxkyz_lo%llim_proc, kxkyz_lo%ulim_proc
              iz = iz_idx(kxkyz_lo, ikxkyz)
              it = it_idx(kxkyz_lo, ikxkyz)
              ikx = ikx_idx(kxkyz_lo, ikxkyz)
              iky = iky_idx(kxkyz_lo, ikxkyz)
              is = is_idx(kxkyz_lo, ikxkyz)
            
              ! Integrate g to get - 2 beta sum_s n_s T_s J1 mu g and store in bpar.
              ! TO DO - Is spec(is)%z supposed to be here?
              call gyro_average_j1(spread(mu, 1, nvpa) * g(:, :, ikxkyz), ikxkyz, g0)
              wgt = -2.0 * beta * spec(is)%z * spec(is)%dens_psi0 * spec(is)%temp_psi0
              call integrate_vmu(g0, iz, tmp)
              bpar(iky, ikx, iz, it) = bpar(iky, ikx, iz, it) + wgt * tmp
          end do
         
          ! Sum the values on all processors and send them to <proc0>.
          call sum_allreduce(bpar)
          deallocate (g0)
      end if

      ! ========================================================================================================================================================= !
      ! ---------- Now that we have the moments of gneo, we need to get the correct fields. This will depend on which fields we are running with. --------------- ! 
      ! ========================================================================================================================================================= !

      ! If we are running electrostatic, we only need to get phi and there is no coupling between fields. Calculate phi by dividing with denominator_fields_neo[iky,ikz,iz].  
      if (fphi > epsilon(0.0) .and. .not. include_apar .and. .not. include_bpar) then
          call calculate_neo_phi(phi, dist, skip_fsa_local)
      end if

      ! If we are running electromagnetic but no bpar, the phi and apar field equations are a 2 x 2 matrix solve. 
      if (fphi > epsilon(0.0) .and. include_apar .and. .not. include_bpar) then
          call calculate_neo_phi_and_apar(phi, apar, dist, skip_fsa_local)
      end if

      ! If we are running electromagnetic fully electromagnetic, the phi apar and bpar field equations are a 3 x 3 matrix.
      if (fphi > epsilon(0.0) .and. include_apar .and. include_bpar) then
          call calculate_neo_phi_apar_and_bpar(phi, apar, bpar, dist, skip_fsa_local)
      end if

   end subroutine advance_fields_fluxtube_using_neo_field_equations_kxkyzlo


! ================================================================================================================================================================================ !
! -------------------------------------------------------------------------- Calculate NEO phi. ---------------------------------------------------------------------------------- !
! ================================================================================================================================================================================ !

    subroutine calculate_neo_phi(phi, dist, skip_fsa)
        ! Parallelisation.
        use mp, only: proc0, mp_abort
        use job_manage, only: time_message
        use multibox, only: mb_calculate_phi
        use timers, only: time_field_solve
      
        ! Arrays.
        use arrays, only: denominator_fields
        use arrays, only: denominator_fields_MBR
        use arrays, only: denominator_fields_h
        use arrays, only: denominator_fields_MBR_h
      
        ! Grids.
        use grids_z, only: nzgrid, ntubes
        use grids_kxky, only: zonal_mode
        use grids_kxky, only: nakx, naky
        use grids_species, only: spec
      
        ! Adiabatic electrons.
        use grids_species, only: has_electron_species
        use grids_species, only: adiabatic_option_switch
        use grids_species, only: adiabatic_option_fieldlineavg
      
        ! Geometry.
        use geometry, only: dl_over_b

        ! HO corrections. 
        use arrays, only: denominator_fields_neo

        implicit none

        ! Arguments.
        character(*), intent(in) :: dist
        complex, dimension(:, :, -nzgrid:, :), intent(in out) :: phi
        logical, optional, intent(in) :: skip_fsa
      
        ! Local variables.
        real, dimension(:, :, :, :), allocatable :: denominator_fields_t_neo
        integer :: ia, it, ikx
        complex :: tmp
        logical :: skip_fsa_local
        logical :: has_elec, adia_elec

        ! ============================================================== !
      
        ! Used for the Dougherty collision operator.
        skip_fsa_local = .false.
        if (present(skip_fsa)) skip_fsa_local = skip_fsa

        ! Assume we only have one field line.
        ia = 1
      
        ! Check if we need to add a Boltzmann response for an adiabatic species.
        has_elec = has_electron_species(spec)
        adia_elec = .not. has_elec .and. (adiabatic_option_switch == adiabatic_option_fieldlineavg)

        ! Start timer.
        if (proc0) call time_message(.false., time_field_solve(:, 4), ' calculate_phi')

        ! ================================================================================================================================================= !
        ! --------------------------------------------------- Using g_neo as the distribution function. --------------------------------------------------- !
        ! ================================================================================================================================================= !
        !                                                                                                                                                   !
        ! If we are using the g_neo distribution, then:                                                                                                     ! 
        !                                                                                                                                                   ! 
        !     phi = sum_s Z_s n_s [ (2B/sqrt(pi)) int dvpa int dmu J_0 * g ]                                                                                ! 
        !     / [ sum_s (Z_s² n_s/T_s) (2B/sqrt(pi)) int dvpa int dmu exp(-v²) (1 - J_0^2) (1 + F_1 - dH_1/dμ|_v∥ 0.5/B_0)]                                 ! 
        !                                                                                                                                                   ! 
        ! denominator_fields_neo[iky,ikz,iz] = [ sum_s (Z_s² n_s/T_s) (2B/sqrt(pi)) int dvpa int dmu exp(-v²) (1 - J_0^2) (1 + F_1 - dH_1/dμ|_v∥ 0.5/B_0)]  ! 
        !                                                                                                                                                   !
        ! To avoid any issues with division by zero we set phi = 0.0 if the denominator is too small. Only thing that makes sense is to set phi = 0.0 if    !
        ! the prefactor for phi in QN is also zero.                                                                                                         !
        !                                                                                                                                                   !
        ! ================================================================================================================================================= !     

        if (dist == 'gneo') then  
            allocate (denominator_fields_t_neo(naky, nakx, -nzgrid:nzgrid, ntubes))
            denominator_fields_t_neo = spread(denominator_fields_neo, 4, ntubes)
            where (denominator_fields_t_neo < epsilon(0.0))
                phi = 0.0
            elsewhere
                phi = phi / denominator_fields_t_neo
            end where
            deallocate (denominator_fields_t_neo)
        else
            ! Abort if <dist> is not recognized.
            if (proc0) write (*, *) 'unknown dist option in calculate_neo_phi. aborting'
            call mp_abort('unknown dist option in calculate_neo_phi. aborting')
            return
        end if

        ! The kx = ky = 0.0 mode is not evolved by stella so make sure this term is set to zero.
        if (any(denominator_fields(1, 1, :) < epsilon(0.))) phi(1, 1, :, :) = 0.0
        if (proc0) call time_message(.false., time_field_solve(:, 4), ' calculate_phi')

        if (adia_elec .and. zonal_mode(1) .and. .not. skip_fsa_local) then
            if (dist == 'gneo') then
                do ikx = 1, nakx
                    do it = 1, ntubes
                        tmp = sum(dl_over_b(ia, :) * phi(1, ikx, :, it))
                        phi(1, ikx, :, it) = phi(1, ikx, :, it) + tmp * denominator_fields_MBR(ikx, :)
                    end do
                end do
            else
                if (proc0) write (*, *) 'unknown dist option in calculate_neo_phi. aborting'
                call mp_abort('unknown dist option in calcuate_neo_phi. aborting')
            end if
        end if
   end subroutine calculate_neo_phi


! ================================================================================================================================================================================ !
! ------------------------------------------------ Calculate phi and apar for electromagnetic simulations when bpar is not included. --------------------------------------------- !
! ================================================================================================================================================================================ !

    subroutine calculate_neo_phi_and_apar(phi, apar, dist, skip_fsa)
        ! Parallelisation.
        use mp, only: proc0, mp_abort
        use job_manage, only: time_message
        use timers, only: time_field_solve

        ! Arrays.
        use arrays, only: denominator_fields_neo, denominator_fields_neo_12
        use arrays, only: denominator_fields_neo_21, denominator_fields_neo_22_g, denominator_fields_neo_22_gbar

        ! Grids.
        use grids_z, only: nzgrid, ntubes
        use grids_kxky, only: nakx, naky

        implicit none

        ! Arguments.
        complex, dimension(:, :, -nzgrid:, :), intent(in out) :: phi, apar
        character(*), intent(in) :: dist
        logical, optional, intent(in) :: skip_fsa

        ! Local variables.
        integer :: ia, it, ikx, iky, iz
        logical :: skip_fsa_local

        ! LAPACK Variables.
        complex(8)        :: A_lapack(2,2)
        complex(8)        :: B_lapack(2,1)
        integer        :: ipiv(2)
        integer        :: info
        external zgesv

        ! ======================================================================================================================================================== ! 
        ! Due to the presence of F_1, all fluctuating fields now couple to one another in the field equations. In the abscence of bpar, this reduces to a 2 x 2    !
        ! matrix problem where the solution provides phi and apar. This is solved with LAPACK. This could be solved directly using Cramer's rule, as is done for   ! 
        ! the coupling of phi and bpar at leading order. However, this would become algebriaclly cumbersome in the fully electromagnetic case where the matrix     !
        ! becomes 3 x 3. The extension of the LAPACK logic to the fully electromagnetic regime is by comparison much easier.                                       !
        ! ======================================================================================================================================================== !

        ! Used for the Dougherty collision operator.
        skip_fsa_local = .false.
        if (present(skip_fsa)) skip_fsa_local = skip_fsa

        ! Assume we only have one field line.
        ia = 1

        if (dist == 'gneo' .or. dist == 'gbarneo') then
            do it = 1, ntubes
                do iz = -nzgrid, nzgrid
                    do ikx = 1, nakx
                        do iky = 1, naky
                            ! Promote real matrix elements to complex numbers to be solved with zgesv. 
                            A_lapack(1,1) = cmplx(denominator_fields_neo(iky,ikx,iz), 0.0)
                            A_lapack(2,1) = cmplx(denominator_fields_neo_21(iky,ikx,iz), 0.0)
                            A_lapack(1,2) = cmplx(denominator_fields_neo_12(iky,ikx,iz), 0.0)
                            if (dist == 'gneo') then
                                A_lapack(2,2) = cmplx(denominator_fields_neo_22_g(iky,ikx,iz), 0.0)
                            else if (dist == 'gbarneo') then
                                A_lapack(2,2) = cmplx(denominator_fields_neo_22_gbar(iky,ikx,iz), 0.0)
                            end if

                            B_lapack(1,1) = phi(iky,ikx,iz,it)
                            B_lapack(2,1) = apar(iky,ikx,iz,it)

                            call zgesv(2, 1, A_lapack, 2, ipiv, B_lapack, 2, info)

                            if (info == 0) then
                                ! Assign solutions to the fields. 
                                phi(iky,ikx,iz,it)  = B_lapack(1,1)
                                apar(iky,ikx,iz,it) = B_lapack(2,1)                                
                            else
                                if (proc0) write(*,*) 'WARNING: ill-conditioned matrix at iky,ikx,iz=', iky, ikx, iz
                                phi(iky,ikx,iz,it)  = cmplx(0.0, 0.0)
                                apar(iky,ikx,iz,it) = cmplx(0.0, 0.0)
                            end if
                        end do
                    end do
                end do
            end do
        else
            if (proc0) write (*, *) 'Unknown dist option in calculate_neo_phi_and_apar. Aborting.'
            call mp_abort('Unknown dist option in calculate_neo_phi_and_apar. Aborting.')
            return
        end if
    end subroutine calculate_neo_phi_and_apar


! ================================================================================================================================================================================ !
! ------------------------------------------------------- Calculate phi, apar and bpar for fully electromagnetic simulations. ---------------------------------------------------- !
! ================================================================================================================================================================================ !

    subroutine calculate_neo_phi_apar_and_bpar(phi, apar, bpar, dist, skip_fsa)
        ! Parallelisation.
        use mp, only: proc0, mp_abort
        use job_manage, only: time_message
        use timers, only: time_field_solve

        ! Arrays.
        use arrays, only: denominator_fields_neo, denominator_fields_neo_12, denominator_fields_neo_13
        use arrays, only: denominator_fields_neo_21, denominator_fields_neo_22_g, denominator_fields_neo_22_gbar, denominator_fields_neo_23
        use arrays, only: denominator_fields_neo_31, denominator_fields_neo_32, denominator_fields_neo_33

        ! Grids.
        use grids_z, only: nzgrid, ntubes
        use grids_kxky, only: nakx, naky

        implicit none

        ! Arguments.
        complex, dimension(:, :, -nzgrid:, :), intent(in out) :: phi, apar, bpar
        character(*), intent(in) :: dist
        logical, optional, intent(in) :: skip_fsa

        ! Local variables.
        integer :: ia, it, ikx, iky, iz
        logical :: skip_fsa_local

        ! LAPACK Variables.
        complex(8)        :: A_lapack(3,3)
        complex(8)        :: B_lapack(3,1)
        integer           :: ipiv(3)
        integer           :: info
        external zgesv

        ! ======================================================================================================================================================== ! 
        ! Due to the presence of F_1, all fluctuating fields now couple to one another in the field equations. In the presence of bpar, this involves a 3 x 3      !
        ! matrix problem where the solution provides phi, apar and bpar. This is solved with LAPACK.                                                               !
        ! ======================================================================================================================================================== !

        ! Used for the Dougherty collision operator.
        skip_fsa_local = .false.
        if (present(skip_fsa)) skip_fsa_local = skip_fsa

        ! Assume we only have one field line.
        ia = 1

        if (dist == 'gneo' .or. dist == 'gbarneo') then
            do it = 1, ntubes
                do iz = -nzgrid, nzgrid
                    do ikx = 1, nakx
                        do iky = 1, naky
                            ! Promote real matrix elements to complex numbers to be solved with zgesv. 
                            A_lapack(1,1) = cmplx(denominator_fields_neo(iky,ikx,iz), 0.0)
                            A_lapack(2,1) = cmplx(denominator_fields_neo_21(iky,ikx,iz), 0.0)
                            A_lapack(3,1) = cmplx(denominator_fields_neo_31(iky,ikx,iz), 0.0)
                            A_lapack(1,2) = cmplx(denominator_fields_neo_12(iky,ikx,iz), 0.0)
                            if (dist == 'gneo') then
                                A_lapack(2,2) = cmplx(denominator_fields_neo_22_g(iky,ikx,iz), 0.0)
                            else if (dist == 'gbarneo') then
                                A_lapack(2,2) = cmplx(denominator_fields_neo_22_gbar(iky,ikx,iz), 0.0)
                            end if
                            A_lapack(3,2) = cmplx(denominator_fields_neo_32(iky,ikx,iz), 0.0)
                            A_lapack(1,3) = cmplx(denominator_fields_neo_13(iky,ikx,iz), 0.0)
                            A_lapack(2,3) = cmplx(denominator_fields_neo_23(iky,ikx,iz), 0.0)
                            A_lapack(3,3) = cmplx(denominator_fields_neo_33(iky,ikx,iz), 0.0)

                            B_lapack(1,1) = phi(iky,ikx,iz,it)
                            B_lapack(2,1) = apar(iky,ikx,iz,it)
                            B_lapack(3,1) = bpar(iky,ikx,iz,it)

                            call zgesv(3, 1, A_lapack, 3, ipiv, B_lapack, 3, info)

                            if (info == 0) then
                                ! Assign solutions to the fields. 
                                phi(iky,ikx,iz,it)  = B_lapack(1,1)
                                apar(iky,ikx,iz,it) = B_lapack(2,1)
                                bpar(iky,ikx,iz,it) = B_lapack(3,1)                                
                            else
                                if (proc0) write(*,*) 'WARNING: ill-conditioned matrix at iky,ikx,iz=', iky, ikx, iz
                                phi(iky,ikx,iz,it)  = cmplx(0.0, 0.0)
                                apar(iky,ikx,iz,it) = cmplx(0.0, 0.0)
                                bpar(iky,ikx,iz,it) = cmplx(0.0, 0.0)
                            end if
                        end do
                    end do
                end do
            end do
        else
            if (proc0) write (*, *) 'Unknown dist option in calculate_neo_phi_and_apar. Aborting.'
            call mp_abort('Unknown dist option in calculate_neo_phi_and_apar. Aborting.')
            return
        end if
    end subroutine calculate_neo_phi_apar_and_bpar


! ================================================================================================================================================================================ !
! -------------------------------------------------- Advances only apar; needed for the implicit advance of the mirror term.  ---------------------------------------------------- !
! ================================================================================================================================================================================ !

    subroutine advance_apar_neo(g, dist, apar)
        ! Parallelisation.
        use mp, only: mp_abort, sum_allreduce, proc0
        use parallelisation_layouts, only: kxkyz_lo
        use parallelisation_layouts, only: iky_idx, ikx_idx, iz_idx, it_idx, is_idx
      
        ! Parameters. 
        use parameters_physics, only: fphi
        use parameters_physics, only: include_apar, include_bpar
        use parameters_physics, only: beta
      
        ! Grids.
        use grids_kxky, only: nakx, naky 
        use grids_species, only: spec
        use grids_z, only: nzgrid, ntubes
        use grids_velocity, only: vpa, mu, nvpa, nmu
        use calculations_velocity_integrals, only: integrate_vmu
      
        ! Calculations.
        use calculations_gyro_averages, only: gyro_average, gyro_average_j1

        ! Arrays. 
        use arrays, only: denominator_fields_neo, denominator_fields_neo_12, denominator_fields_neo_13
        use arrays, only: denominator_fields_neo_21, denominator_fields_neo_22_g, denominator_fields_neo_23
        use arrays, only: denominator_fields_neo_31, denominator_fields_neo_32, denominator_fields_neo_33

        implicit none

        ! Arguments.
        character(*), intent(in) :: dist
        complex, dimension(:, :, kxkyz_lo%llim_proc:), intent(in) :: g
        complex, dimension(:, :, -nzgrid:, :), intent(out) :: apar

        ! Local variables.
        integer :: ikxkyz, iky, ikx, iz, it, is
        integer :: ia
        real :: wgt
        complex :: tmp
        complex, dimension(:, :), allocatable :: scratch
        complex, dimension(:, :, :, :), allocatable :: phi, bpar

        ia = 1    

        ! ==================================================================== !
        
        ! Allocate temporary arrays for the phi and bpar.
        allocate(phi(naky, nakx, -nzgrid:nzgrid, ntubes))
        allocate(bpar(naky, nakx, -nzgrid:nzgrid, ntubes))      

        ! Initialise the phi field.
        phi = 0.

        ! Calculate the phi if it is included.
        if (fphi > epsilon(0.0)) then        
            ! Allocate temporary array.
            allocate (scratch(nvpa, nmu))
          
            ! Iterate over the (kx,ky,z,mu,vpa,s) points.
            do ikxkyz = kxkyz_lo%llim_proc, kxkyz_lo%ulim_proc
                iz = iz_idx(kxkyz_lo, ikxkyz)
                it = it_idx(kxkyz_lo, ikxkyz)
                ikx = ikx_idx(kxkyz_lo, ikxkyz)
                iky = iky_idx(kxkyz_lo, ikxkyz)
                is = is_idx(kxkyz_lo, ikxkyz)
            
                ! Integrate g to get sum_s Z_s n_s J0 g and store in phi.
                call gyro_average(g(:, :, ikxkyz), ikxkyz, scratch)
                wgt = spec(is)%z * spec(is)%dens_psi0
                call integrate_vmu(scratch, iz, tmp)
                phi(iky, ikx, iz, it) = phi(iky, ikx, iz, it) + wgt * tmp
            end do
         
            ! Sum the values on all processors and send them to <proc0>
            call sum_allreduce(phi)    

            ! Deallocate temporary arrays.
            deallocate (scratch)
        end if

        ! Initialise the apar field.
        apar = 0.
      
        if (include_apar) then
            ! Allocate temporary array.
            allocate (scratch(nvpa, nmu))
         
            ! Iterate over the (kx,ky,z,mu,vpa,s) points.
            do ikxkyz = kxkyz_lo%llim_proc, kxkyz_lo%ulim_proc
                iz = iz_idx(kxkyz_lo, ikxkyz)
                it = it_idx(kxkyz_lo, ikxkyz)
                ikx = ikx_idx(kxkyz_lo, ikxkyz)
                iky = iky_idx(kxkyz_lo, ikxkyz)
                is = is_idx(kxkyz_lo, ikxkyz)
            
                call gyro_average(spread(vpa, 2, nmu) * g(:, :, ikxkyz), ikxkyz, scratch)
                wgt = beta * spec(is)%z * spec(is)%dens_psi0 * spec(is)%stm_psi0
                call integrate_vmu(scratch, iz, tmp)
                apar(iky, ikx, iz, it) = apar(iky, ikx, iz, it) + tmp * wgt
            end do
         
            ! Apar for different species may be spread over processors at this point, so broadcast to all procs and sum over species.
            call sum_allreduce(apar)
            
            deallocate (scratch)
        end if

        bpar = 0.
        
        ! Calculate the bpar moment. 
        if (include_bpar) then
            ! Allocate temporary array.
            allocate (scratch(nvpa, nmu))

            ! Iterate over the (kx,ky,z,mu,vpa,s) points.
            do ikxkyz = kxkyz_lo%llim_proc, kxkyz_lo%ulim_proc
                iz = iz_idx(kxkyz_lo, ikxkyz)
                it = it_idx(kxkyz_lo, ikxkyz)
                ikx = ikx_idx(kxkyz_lo, ikxkyz)
                iky = iky_idx(kxkyz_lo, ikxkyz)
                is = is_idx(kxkyz_lo, ikxkyz)

                ! Integrate g to get - 2 beta sum_s n_s T_s J1 mu g and store in bpar.
                call gyro_average_j1(spread(mu, 1, nvpa) * g(:, :, ikxkyz), ikxkyz, scratch)
                wgt = -2.0 * beta * spec(is)%dens_psi0 * spec(is)%temp_psi0
                call integrate_vmu(scratch, iz, tmp)
                bpar(iky, ikx, iz, it) = bpar(iky, ikx, iz, it) + wgt * tmp
            end do

            ! Sum the values on all processors and send them to <proc0>.
            call sum_allreduce(bpar)
            deallocate (scratch)
        end if

        ! With the moments, calculate the fields, but return only apar. 
        call get_apar_neo(phi, apar, bpar, dist)

        ! Deallocate temporary arrays for phi and bpar.
        if (allocated(phi)) deallocate(phi)
        if (allocated(bpar)) deallocate(bpar)      

    end subroutine advance_apar_neo


! ======================================================================================================================================================================================= !
! --------------------------------------------------------------- Gets the correct apar given the fields included in the simuation ------------------------------------------------------ !
! ======================================================================================================================================================================================= !

    subroutine get_apar_neo(phi, apar, bpar, dist)
        ! Parallelisation.
        use mp, only: proc0, mp_abort
      
        ! Arrays.
        use arrays, only: denominator_fields_neo, denominator_fields_neo_12, denominator_fields_neo_13
        use arrays, only: denominator_fields_neo_21, denominator_fields_neo_22_g, denominator_fields_neo_23
        use arrays, only: denominator_fields_neo_31, denominator_fields_neo_32, denominator_fields_neo_33

        ! Parameters.
        use parameters_physics, only: include_apar, include_bpar
      
        ! Grids
        use grids_z, only: nzgrid, ntubes
        use grids_kxky, only: nakx, naky 
      
        implicit none

        ! Arguments.
        complex, dimension(:, :, -nzgrid:, :), intent(in) :: phi, bpar
        complex, dimension(:, :, -nzgrid:, :), intent(in out) :: apar
        character(*), intent(in) :: dist

        ! Local variables.
        integer :: ia
        integer :: ikxkyz, iky, ikx, iz, it, is      

        ! LAPACK Variables.
        complex(8)     :: A_lapack(2,2)
        complex(8)     :: B_lapack(2,1)
        complex(8)     :: C_lapack(3,3)
        complex(8)     :: D_lapack(3,1)
        integer        :: ipiv(2), jpiv(3)
        integer        :: info
        external zgesv

        ! Assume we only have one field line
        ia = 1

        ! ================================================================================================================= !
      
        ! We have the sources that we need to get the associated field matrix for apar.
        ! Note that the implicit mirror advance only uses gneo and so there is no gbarneo dist option here.
        ! If we are only including apar, only solve the 2x2 matrix.  
        if (include_apar .and. .not. include_bpar) then 
            if (dist == 'gneo') then
                do it = 1, ntubes
                    do iz = -nzgrid, nzgrid
                        do ikx = 1, nakx
                            do iky = 1, naky
                                A_lapack(1,1) = cmplx(denominator_fields_neo(iky,ikx,iz), 0.0)
                                A_lapack(2,1) = cmplx(denominator_fields_neo_21(iky,ikx,iz), 0.0)
                                A_lapack(1,2) = cmplx(denominator_fields_neo_12(iky,ikx,iz), 0.0)
                                A_lapack(2,2) = cmplx(denominator_fields_neo_22_g(iky,ikx,iz), 0.0)

                                B_lapack(1,1) = phi(iky,ikx,iz,it)
                                B_lapack(2,1) = apar(iky,ikx,iz,it)
     
                                call zgesv(2, 1, A_lapack, 2, ipiv, B_lapack, 2, info)

                                if (info == 0) then 
                                    apar(iky,ikx,iz,it) = B_lapack(2,1)
                                else
                                    if (proc0) write(*,*) 'WARNING: ill-conditioned matrix in get_apar_neo at iky, ikx, iz, it =', iky, ikx, iz, it
                                    apar(iky,ikx,iz,it) = cmplx(0.0, 0.0)
                                end if
                            end do
                        end do
                    end do
                end do
            else  
                if (proc0) write (*, *) 'Unknown dist option in get_apar_neo. Aborting.'
                call mp_abort('Unknown dist option in get_apar_neo. Aborting.')
                return      
            end if     
        end if

        ! If we are including apar and bpar, solve the full 3x3 matrix. 
        if (include_apar .and. include_bpar) then
            if (dist == 'gneo') then
                do it = 1, ntubes
                    do iz = -nzgrid, nzgrid
                        do ikx = 1, nakx
                            do iky = 1, naky
                                C_lapack(1,1) = cmplx(denominator_fields_neo(iky,ikx,iz), 0.0)
                                C_lapack(2,1) = cmplx(denominator_fields_neo_21(iky,ikx,iz), 0.0)
                                C_lapack(3,1) = cmplx(denominator_fields_neo_31(iky,ikx,iz), 0.0)
                                C_lapack(1,2) = cmplx(denominator_fields_neo_12(iky,ikx,iz), 0.0)
                                C_lapack(2,2) = cmplx(denominator_fields_neo_22_g(iky,ikx,iz), 0.0)
                                C_lapack(3,2) = cmplx(denominator_fields_neo_32(iky,ikx,iz), 0.0)
                                C_lapack(1,3) = cmplx(denominator_fields_neo_13(iky,ikx,iz), 0.0)
                                C_lapack(2,3) = cmplx(denominator_fields_neo_23(iky,ikx,iz), 0.0)
                                C_lapack(3,3) = cmplx(denominator_fields_neo_33(iky,ikx,iz), 0.0)

                                D_lapack(1,1) = phi(iky,ikx,iz,it)
                                D_lapack(2,1) = apar(iky,ikx,iz,it)
                                D_lapack(3,1) = bpar(iky,ikx,iz,it)

                                call zgesv(3, 1, C_lapack, 3, jpiv, D_lapack, 3, info)

                                if (info == 0) then
                                    apar(iky,ikx,iz,it) = D_lapack(2,1)
                                else
                                    if (proc0) write(*,*) 'WARNING: ill-conditioned matrix in get_apar_neo at iky, ikx, iz, it =', iky, ikx, iz, it
                                    apar(iky,ikx,iz,it) = cmplx(0.0, 0.0)
                                end if
                            end do
                        end do
                    end do
                end do
            else
                if (proc0) write (*, *) 'Unknown dist option in get_apar_neo. Aborting.'
                call mp_abort('Unknown dist option in get_apar_neo. Aborting.')
                return
            end if
        end if
   end subroutine get_apar_neo


! ================================================================================================================================================================================ !
! --------------------------------------------  Provides the extended denominator for phi when running electrostatic HO simulations ---------------------------------------------- !
! ================================================================================================================================================================================ !

   subroutine init_neo_electrostatic_fields
      ! Parallelisation.
      use mp, only: sum_allreduce
      use parallelisation_layouts, only: kxkyz_lo, iz_idx, it_idx, ikx_idx, iky_idx, is_idx
      
      ! Arrays.
      use arrays, only: denominator_fields_neo, efac
      use arrays_gyro_averages, only: aj0v

      ! Parameters.
      use parameters_physics, only: fphi
      
      ! Grids.
      use grids_velocity, only: nvpa, nmu
      use grids_species, only: spec
      use grids_kxky, only: zonal_mode, akx
      use grids_kxky, only: nakx
      use grids_velocity, only: maxwell_vpa, maxwell_mu, maxwell_fac
      use grids_z, only: nzgrid

      ! Adiabatic electrons.
      use grids_species, only: has_electron_species
      use grids_species, only: ion_species
      use grids_species, only: tite, nine
      use grids_species, only: adiabatic_option_switch
      use grids_species, only: adiabatic_option_fieldlineavg

      ! Calculations.
      use calculations_velocity_integrals, only: integrate_vmu
      
      ! Geometry.
      use geometry, only: dl_over_b, bmag

      ! NEO data. 
      use neoclassical_terms_neo, only: neo_mu_fac_global

      implicit none

      ! Local variables.
      integer :: ikxkyz, iz, it, ikx, iky, ia, is
      real :: tmp, wgt
      real, dimension(:, :), allocatable :: g0

      ! Assume we only have a single field line.
      ia = 1

      ! Allocate denominator_fields_neo. 
      call allocate_neo_electrostatic_fields

      ! Calculate the denominators needed for electrostatic neoclassical simulations.
      if (fphi > epsilon(0.0)) then
      
         ! Allocate temporary arrays.
         allocate (g0(nvpa, nmu))

         ! ======================================================================================================================================================== ! 
         ! When we are evolving NEO's higher order corrections, we use the distribution function g_neo. The corresponding electrostatic QN condition for phi is:    ! 
         !                                                                                                                                                          ! 
         !     phi = sum_s Z_s n_s [ (2B/sqrt(pi)) int dvpa int dmu J_0 * g_neo ]                                                                                       ! 
         !     / [ sum_s (Z_s² n_s/T_s) (2B/sqrt(pi)) int dvpa int dmu exp(-v²) (1 - J_0^2) (1 + F_1 - dH_1/dμ|_v∥ 0.5/B_0)]                                        ! 
         !                                                                                                                                                          ! 
         !     denominator_fields_neo[iky,ikz,iz] = [ sum_s (Z_s² n_s/T_s) (2B/sqrt(pi)) int dvpa int dmu exp(-v²) (1 - J_0^2) (1 + F_1 - dH_1/dμ|_v∥ 0.5/B_0)]     ! 
         !                                                                                                                                                          ! 
         ! ======================================================================================================================================================== !
         
         do ikxkyz = kxkyz_lo%llim_proc, kxkyz_lo%ulim_proc
            ! <denominator_fields_neo> does not depend on flux tube index, so only compute for one flux tube index.
            it = it_idx(kxkyz_lo, ikxkyz)
            if (it /= 1) cycle
            iky = iky_idx(kxkyz_lo, ikxkyz)
            ikx = ikx_idx(kxkyz_lo, ikxkyz)
            iz = iz_idx(kxkyz_lo, ikxkyz)
            is = is_idx(kxkyz_lo, ikxkyz)

            ! Calculate (1 - J_0(a_k)²) exp(v²) for each (kx,ky,z).
            g0 = spread((1.0 - aj0v(:, ikxkyz)**2), 1, nvpa) * spread(maxwell_vpa(:, is), 2, nmu) * spread(maxwell_mu(ia, iz, :, is), 1, nvpa) * maxwell_fac(is)

            ! Multiply this by (1 + F_1 - dH_1/dμ|_v∥ 0.5/B_0) for each (kx,ky,z). 
            g0 = g0 * ( 1.0 - 0.5 * neo_mu_fac_global(iz, :, :, is, 1) / bmag(ia, iz) )
            
            ! Calculate denominator_fields_neo[iky,ikz,iz].
            wgt = spec(is)%z * spec(is)%z * spec(is)%dens_psi0 / spec(is)%temp
            call integrate_vmu(g0, iz, tmp)
            denominator_fields_neo(iky, ikx, iz) = denominator_fields_neo(iky, ikx, iz) + tmp * wgt
         end do
         
         ! Sum the values on all processors and send them to <proc0>.
         call sum_allreduce(denominator_fields_neo)
         
         ! Avoid divide by zero when kx=ky=0; We do not evolve this mode, so the value is irrelevant.
         if (zonal_mode(1) .and. akx(1) < epsilon(0.) .and. has_electron_species(spec)) then
            denominator_fields_neo(1, 1, :) = 0.0
         end if

         ! ======================================================================================================================================================== ! 
         !                                                                                                                                                          ! 
         ! When using adiabatic electrons, denominator_fields_neo acquires a factor associated with the adiabatic response:                                         !
         !                                                                                                                                                          ! 
         ! denominator_fields_neo[iky,ikz,iz] = denominator_fields_neo[iky,ikz,iz] +2B/sqrt(pi) int dvpa int dmu exp(-v²) (F_{1,e} - Z_e * ϕ^1_0)                   !
         !                                                                                                                                                          !
         ! This factor is the lowest order moment of the electron F_1 distribution and is computed in neoclassical_terms_neo.f90 as "neo_dens".                     !
         ! Here is it just added to the denominator.                                                                                                                !
         !                                                                                                                                                          !
         ! ======================================================================================================================================================== !

         if (.not. has_electron_species(spec)) then
             denominator_fields_neo = denominator_fields_neo + efac
         end if

         ! Deallocate temporary arrays.
         if (allocated(g0)) deallocate (g0)
      end if
   end subroutine init_neo_electrostatic_fields


! ================================================================================================================================================================================ !
! --------------------------------------- Provides the matrix elements needed for phi, apar and bpar when running electromagnetic HO simulations. -------------------------------- !
! ================================================================================================================================================================================ !

    subroutine init_neo_electromagnetic_fields
        ! Parallelisation.
        use mp, only: sum_allreduce
        use parallelisation_layouts, only: kxkyz_lo
        use parallelisation_layouts, only: iz_idx, it_idx, ikx_idx, iky_idx, is_idx

        ! Arrays.
        use arrays, only: kperp2
        use arrays, only: denominator_fields_neo_12, denominator_fields_neo_13
        use arrays, only: denominator_fields_neo_21, denominator_fields_neo_22_g, denominator_fields_neo_22_gbar, denominator_fields_neo_23
        use arrays, only: denominator_fields_neo_31, denominator_fields_neo_32, denominator_fields_neo_33

        ! Parameters.
        use parameters_physics, only: include_apar, include_bpar, beta
        use parameters_physics, only: fphi
      
        ! Geometry.
        use geometry, only: bmag

        ! Grids.
        use grids_kxky, only : nakx, naky
        use grids_z, only: nzgrid
        use grids_species, only: spec
        use grids_velocity, only: vpa, mu, nvpa, nmu
        use grids_velocity, only: maxwell_vpa, maxwell_mu, maxwell_fac
      
        ! Calculations.
        use calculations_velocity_integrals, only: integrate_vmu
        use arrays_gyro_averages, only: aj0v, aj1v

        ! Neoclassical data.
        use neoclassical_terms_neo, only: neo_mu_fac_global, neo_vpa_fac_global

        implicit none

        ! Local variables.
        integer :: ikxkyz, iz, it, ikx, iky, is, ia, iv, imu
        real :: tmp, wgt, denom_tmp
        real, dimension(:, :), allocatable :: g0
        
        if (.not. (include_apar .or. include_bpar)) return

        ! Allocate the arrays needed for higher order electromagnetic field solve calculations.
        call allocate_neo_electromagnetic_fields

        ia = 1 
    
        ! ======================================================================================================================================================== ! 
        ! denominator_fields_neo_12 is the apar contribution to the QN condition. This is given by:                                                                ! 
        !                                                                                                                                                          ! 
        ! denominator_fields_neo_12[iky,ikz,iz] = { sum_s (Z_s² v_{th,s} n_s/T_s) (2B/sqrt(pi)) int dvpa int dmu exp(-v²)                                          ! 
        ! [vpa * (1 - J_0^2) (dH_1/dμ|_v∥ - 2B_0 * F_1) / B_0  + (2 vpa * F_1 - dH_1/dv∥|_μ)]}                                                                             !  
        !                                                                                                                                                          ! 	 
        ! ======================================================================================================================================================== !

        if (include_apar) then 
            allocate (g0(nvpa, nmu))
        
            do ikxkyz = kxkyz_lo%llim_proc, kxkyz_lo%ulim_proc
                it = it_idx(kxkyz_lo, ikxkyz)
                if (it /= 1) cycle
                iky = iky_idx(kxkyz_lo, ikxkyz)
                ikx = ikx_idx(kxkyz_lo, ikxkyz)
                iz = iz_idx(kxkyz_lo, ikxkyz)
                is = is_idx(kxkyz_lo, ikxkyz)

                ! Calculate the Maxwellian factor.
                g0 = spread(maxwell_vpa(:, is), 2, nmu) * maxwell_fac(is) * spread(maxwell_mu(ia, iz, :, is), 1, nvpa)

                ! Multiply by the neoclassical factor.
                g0 = g0 * ( ( spread((1.0 - aj0v(:, ikxkyz)**2), 1, nvpa) * spread(vpa, 2, nmu) * neo_mu_fac_global(iz, :, :, is, 1) / bmag(ia, iz) ) &
                - neo_vpa_fac_global(iz, :, :, is, 1) )

                ! Calculate denominator_fields_neo_12[iky,ikz,iz].
                wgt = spec(is)%z * spec(is)%z * spec(is)%dens_psi0 * spec(is)%stm / spec(is)%temp 
                call integrate_vmu(g0, iz, tmp)
                denominator_fields_neo_12(iky, ikx, iz) = denominator_fields_neo_12(iky, ikx, iz) + tmp * wgt
            end do

            ! Sum the values on all processors and send them to <proc0>.
            call sum_allreduce(denominator_fields_neo_12)

            ! ======================================================================================================================================================== ! 
            ! denominator_fields_neo_21 is the phi contribution to parallel Amperes law. This is given by:                                                             ! 
            !                                                                                                                                                          ! 
            ! denominator_fields_neo_21[iky,ikz,iz] = β sum_s (Z_s² v_{th,s} n_s/T_s) (2B/sqrt(pi)) int dvpa int dmu exp(-v²) vpa                                      ! 
            ! [(1 - J_0^2) ( F_1 - 0.5 * dH_1/dμ|_v∥/B_0 )                                                                                                             !  
            !                                                                                                                                                          ! 
            ! ======================================================================================================================================================== !

            do ikxkyz = kxkyz_lo%llim_proc, kxkyz_lo%ulim_proc
                it = it_idx(kxkyz_lo, ikxkyz)
                if (it /= 1) cycle
                iky = iky_idx(kxkyz_lo, ikxkyz)
                ikx = ikx_idx(kxkyz_lo, ikxkyz)
                iz = iz_idx(kxkyz_lo, ikxkyz)
                is = is_idx(kxkyz_lo, ikxkyz)                 

                ! Calculate Maxwellian x vpa factor.
                g0 = spread(maxwell_vpa(:, is) * vpa, 2, nmu) * maxwell_fac(is) * spread(maxwell_mu(ia, iz, :, is), 1, nvpa) 

                g0 = - g0 * 0.5 * spread((1.0 - aj0v(:, ikxkyz)**2), 1, nvpa) *  neo_mu_fac_global(iz, :, :, is, 1) / bmag(ia, iz)
                
                ! Calculate denominator_fields_neo_21[iky,ikz,iz].
                wgt = beta * spec(is)%z * spec(is)%z * spec(is)%dens_psi0 * spec(is)%stm / spec(is)%temp
                call integrate_vmu(g0, iz, tmp)
                denominator_fields_neo_21(iky, ikx, iz) = denominator_fields_neo_21(iky, ikx, iz) + tmp * wgt
            end do

            ! Sum the values on all processors and send them to <proc0>.
            call sum_allreduce(denominator_fields_neo_21)

            ! ======================================================================================================================================================== ! 
            ! denominator_fields_neo_22_g is the apar contribution to parallel Amperes law when gneo is being used for the distribution function. This is given by:    ! 
            !                                                                                                                                                          !
            ! denominator_fields_neo_22_g[iky,ikz,iz] = kperp2 + β sum_s (Z_s² n_s n_s/m_s) (2B/sqrt(pi)) int dvpa int dmu exp(-v²) vpa                                ! 
            ! { vpa * (1 - J_0^2) * ( dH_1/dμ|_v∥ - 2B_0 F_1 ) / B_0 + ( 2v∥ * F_1 -  dH_1/dv∥|_μ ) }                                                                  !   
            !                                                                                                                                                          !     
            ! ======================================================================================================================================================== !
         
            do ikxkyz = kxkyz_lo%llim_proc, kxkyz_lo%ulim_proc
                it = it_idx(kxkyz_lo, ikxkyz)
                if (it /= 1) cycle
                iky = iky_idx(kxkyz_lo, ikxkyz)
                ikx = ikx_idx(kxkyz_lo, ikxkyz)
                iz = iz_idx(kxkyz_lo, ikxkyz)
                is = is_idx(kxkyz_lo, ikxkyz)

                g0 = spread(maxwell_vpa(:, is)*vpa, 2, nmu) * maxwell_fac(is) * spread(maxwell_mu(ia, iz, :, is), 1, nvpa) 

                g0 = g0 * ( ( spread((1.0 - aj0v(:, ikxkyz)**2), 1, nvpa) * spread(vpa, 2, nmu) * neo_mu_fac_global(iz, :, :, is, 1) / bmag(ia, iz) ) &
                - neo_vpa_fac_global(iz, :, :, is, 1)  )

                ! Calculate denominator_fields_neo_22_g[iky,ikz,iz].
                wgt = beta * spec(is)%z * spec(is)%z * spec(is)%dens / spec(is)%mass
                call integrate_vmu(g0, iz, tmp)
                denominator_fields_neo_22_g(iky, ikx, iz) = denominator_fields_neo_22_g(iky, ikx, iz) + tmp * wgt
            end do

            ! Sum the values on all processors and send them to <proc0>.
            call sum_allreduce(denominator_fields_neo_22_g)

            ! Add the kperp2 factor.
            denominator_fields_neo_22_g = denominator_fields_neo_22_g + kperp2(:, :, ia, :)

            ! ======================================================================================================================================================== ! 
            ! denominator_fields_neo_22_gbar is the apar contribution to parallel Amperes law when gbar_neo is being used for the distribution function. This is       !
            ! given by:                                                                                                                                                ! 
            !                                                                                                                                                          !
            ! denominator_fields_neo_22_gbar[iky,ikz,iz] = kperp2 + β sum_s (Z_s² n_s n_s/m_s) (2B/sqrt(pi)) int dvpa int dmu exp(-v²) vpa                             ! 
            ! { vpa * (1 - J_0^2) * ( dH_1/dμ|_v∥ - 2B_0 F_1 ) / B_0 + ( 2v∥ * F_1 -  dH_1/dv∥|_μ ) }                                                                  !   
            !                                                                                                                                                          !     
            ! ======================================================================================================================================================== !

            do ikxkyz = kxkyz_lo%llim_proc, kxkyz_lo%ulim_proc
                it = it_idx(kxkyz_lo, ikxkyz)
                if (it /= 1) cycle
                iky = iky_idx(kxkyz_lo, ikxkyz)
                ikx = ikx_idx(kxkyz_lo, ikxkyz)
                iz = iz_idx(kxkyz_lo, ikxkyz)
                is = is_idx(kxkyz_lo, ikxkyz)

                g0 = spread(maxwell_vpa(:, is)*vpa, 2, nmu) * maxwell_fac(is) * spread(maxwell_mu(ia, iz, :, is), 1, nvpa) 

                g0 = g0 * ( ( spread((1.0 - aj0v(:, ikxkyz)**2), 1, nvpa) * spread(vpa(:), 2, nmu) * neo_mu_fac_global(iz, :, :, is, 1) / bmag(ia, iz) ) &
                - neo_vpa_fac_global(iz, :, :, is, 1) + 2.0 * spread(vpa, 2, nmu) * spread(aj0v(:, ikxkyz)**2, 1, nvpa) )

                ! Calculate denominator_fields_neo_22_gbar[iky,ikz,iz].            
                wgt = beta * spec(is)%z * spec(is)%z * spec(is)%dens / spec(is)%mass
               
                call integrate_vmu(g0, iz, tmp)
                denominator_fields_neo_22_gbar(iky, ikx, iz) = denominator_fields_neo_22_gbar(iky, ikx, iz) + tmp * wgt
            end do

            ! Sum the values on all processors and send them to <proc0>.
            call sum_allreduce(denominator_fields_neo_22_gbar)

            ! Add the kperp2 factor.
            denominator_fields_neo_22_gbar = denominator_fields_neo_22_gbar + kperp2(:, :, ia, :)

            ! Deallocate temporary array. 
            deallocate (g0)
        end if

        if (include_bpar) then
            ! Allocate temporary array.
            allocate (g0(nvpa, nmu))

            ! ======================================================================================================================================================== ! 
            ! denominator_fields_neo_13 is the bpar contribution to Quasineutrality. This is given by:                                                                 ! 
            !                                                                                                                                                          ! 
            ! denominator_fields_neo_13[iky,ikz,iz] = 4 sum_s Z_s n_s (2B/sqrt(pi)) int dvpa int dmu exp(-v²) mu                                                       !
            ! * J_0 * (J_1 / a_k) * ( 0.5 * dH_1/dμ|_v∥ / B_0 - F_1 - 1 ) / B_0                                                                                        !  
            !                                                                                                                                                          !     
            ! ======================================================================================================================================================== !

            do ikxkyz = kxkyz_lo%llim_proc, kxkyz_lo%ulim_proc
                it = it_idx(kxkyz_lo, ikxkyz)
                if (it /= 1) cycle
                iky = iky_idx(kxkyz_lo, ikxkyz)
                ikx = ikx_idx(kxkyz_lo, ikxkyz)
                iz = iz_idx(kxkyz_lo, ikxkyz)
                is = is_idx(kxkyz_lo, ikxkyz)

                g0 = spread((mu(:) * aj0v(:, ikxkyz) * aj1v(:, ikxkyz)), 1, nvpa) * spread(maxwell_vpa(:, is), 2, nmu) * spread(maxwell_mu(ia, iz, :, is), 1, nvpa) &
                * maxwell_fac(is) * ( ( 0.5 * neo_mu_fac_global(iz, :, :, is, 1) / bmag(ia, iz) ) - 1.0 )
            
                wgt = 4.0 * spec(is)%z * spec(is)%dens_psi0
                call integrate_vmu(g0, iz, tmp)
                denominator_fields_neo_13(iky, ikx, iz) = denominator_fields_neo_13(iky, ikx, iz) + tmp * wgt
            end do

            call sum_allreduce(denominator_fields_neo_13)

            ! ======================================================================================================================================================== ! 
            ! denominator_fields_neo_23 is the bpar contribution to Parallel Ampere's law. This is given by:                                                           ! 
            !                                                                                                                                                          ! 
            ! denominator_fields_neo_23[iky,ikz,iz] = 4β sum_s Z_s n_s v_{th,s} (2B/sqrt(pi)) int dvpa int dmu exp(-v²) mu vpa                                         !
            ! * J_0 * (J_1 / a_k) * ( 0.5 * dH_1/dμ|_v∥ / B_0 - F_1 ) / B_0                                                                                        !  
            !                                                                                                                                                          !     
            ! ======================================================================================================================================================== !

            do ikxkyz = kxkyz_lo%llim_proc, kxkyz_lo%ulim_proc
                it = it_idx(kxkyz_lo, ikxkyz)
                if (it /= 1) cycle
                iky = iky_idx(kxkyz_lo, ikxkyz)
                ikx = ikx_idx(kxkyz_lo, ikxkyz)
                iz = iz_idx(kxkyz_lo, ikxkyz)
                is = is_idx(kxkyz_lo, ikxkyz)

                g0 = spread((mu(:) * aj0v(:, ikxkyz) * aj1v(:, ikxkyz)), 1, nvpa) * spread(vpa(:) * maxwell_vpa(:, is), 2, nmu) * spread(maxwell_mu(ia, iz, :, is), 1, nvpa) &
                * maxwell_fac(is) * neo_mu_fac_global(iz, :, :, is, 1) / bmag(ia, iz)

                wgt = 2.0 * beta *spec(is)%z * spec(is)%dens_psi0 * spec(is)%stm
                call integrate_vmu(g0, iz, tmp)
                denominator_fields_neo_23(iky, ikx, iz) = denominator_fields_neo_23(iky, ikx, iz) + tmp * wgt
            end do

            call sum_allreduce(denominator_fields_neo_23)

            ! ======================================================================================================================================================== ! 
            ! denominator_fields_neo_31 is the phi contribution to Perpendicular Ampere's law. This is given by:                                                       ! 
            !                                                                                                                                                          ! 
            ! denominator_fields_neo_31[iky,ikz,iz] = 2β sum_s Z_s n_s (2B/sqrt(pi)) int dvpa int dmu exp(-v²) mu * { J_0 * (J_1 / a_k)                                 !
            ! + [ 1 - J_0 * (J_1 / a_k) ] * ( 0.5 * dH_1/dμ|_v∥ / B_0 - F_1) }                                                                                         !  
            !                                                                                                                                                          !     
            ! ======================================================================================================================================================== !

            do ikxkyz = kxkyz_lo%llim_proc, kxkyz_lo%ulim_proc
                it = it_idx(kxkyz_lo, ikxkyz)
                if (it /= 1) cycle
                iky = iky_idx(kxkyz_lo, ikxkyz)
                ikx = ikx_idx(kxkyz_lo, ikxkyz)
                iz = iz_idx(kxkyz_lo, ikxkyz)
                is = is_idx(kxkyz_lo, ikxkyz)

                g0 = spread((mu(:) * maxwell_mu(ia, iz, :, is)), 1, nvpa) * spread(maxwell_vpa(:, is), 2, nmu) * maxwell_fac(is) &
                * ( spread((aj0v(:, ikxkyz) * aj1v(:, ikxkyz)), 1, nvpa) + 0.5 * ( 1.0 - spread((aj0v(:, ikxkyz) * aj1v(:, ikxkyz)), 1, nvpa) ) &
                * neo_mu_fac_global(iz, :, :, is, 1) / bmag(ia, iz) )

                wgt = 2.0 * beta * spec(is)%z * spec(is)%dens_psi0
                call integrate_vmu(g0, iz, tmp)
                denominator_fields_neo_31(iky, ikx, iz) = denominator_fields_neo_31(iky, ikx, iz) + tmp * wgt
            end do

            call sum_allreduce(denominator_fields_neo_31)

            ! ======================================================================================================================================================== ! 
            ! denominator_fields_neo_32 is the apar contribution to Perpendicular Ampere's law. This is given by:                                                      ! 
            !                                                                                                                                                          ! 
            ! denominator_fields_neo_32[iky,ikz,iz] = 2β sum_s Z_s n_s v_{th,s} (2B/sqrt(pi)) int dvpa int dmu exp(-v²) mu * { ( dH_1/dv∥|_μ - 2v∥ * F_1 )             !
            ! + v∥ * [ 1 - J_0 * (J_1 / a_k) ] * ( 2B_0 * F_1 - dH_1/dμ|_v∥) }                                                                                         !   
            !                                                                                                                                                          !     
            ! ======================================================================================================================================================== !

            do ikxkyz = kxkyz_lo%llim_proc, kxkyz_lo%ulim_proc
                it = it_idx(kxkyz_lo, ikxkyz)
                if (it /= 1) cycle
                iky = iky_idx(kxkyz_lo, ikxkyz)
                ikx = ikx_idx(kxkyz_lo, ikxkyz)
                iz = iz_idx(kxkyz_lo, ikxkyz)
                is = is_idx(kxkyz_lo, ikxkyz)

                g0 = spread((mu(:) * maxwell_mu(ia, iz, :, is)), 1, nvpa) * spread(maxwell_vpa(:, is), 2, nmu) * maxwell_fac(is) &
                * ( neo_vpa_fac_global(iz, :, :, is, 1) & 
                - spread(vpa(:), 2, nmu) * neo_mu_fac_global(iz, :, :, is, 1) * ( 1.0 - spread((aj0v(:, ikxkyz) * aj1v(:, ikxkyz)), 1, nvpa) ) / bmag(ia, iz) ) 

                wgt = 2.0 * beta * spec(is)%z * spec(is)%dens_psi0 * spec(is)%stm
                call integrate_vmu(g0, iz, tmp)
                denominator_fields_neo_32(iky, ikx, iz) = denominator_fields_neo_32(iky, ikx, iz) + tmp * wgt
            end do

            call sum_allreduce(denominator_fields_neo_32)

            ! ======================================================================================================================================================== ! 
            ! denominator_fields_neo_33 is the bpar contribution to Quasineutrality. This is given by:                                                                 ! 
            !                                                                                                                                                          ! 
            ! denominator_fields_neo_33[iky,ikz,iz] = 1 + 8β sum_s n_s T_s 2B/sqrt(pi) int dvpa int dmu exp(-v²) * mu² * (J_1² / a_k²)                                 ! 
            ! ( 1 + F_1 - 0.5 * dH_1/dμ|_v∥ / B_0 )                                                                                                                       !
            !                                                                                                                                                          !     
            ! ======================================================================================================================================================== !

            do ikxkyz = kxkyz_lo%llim_proc, kxkyz_lo%ulim_proc
                it = it_idx(kxkyz_lo, ikxkyz)
                if (it /= 1) cycle
                iky = iky_idx(kxkyz_lo, ikxkyz)
                ikx = ikx_idx(kxkyz_lo, ikxkyz)
                iz = iz_idx(kxkyz_lo, ikxkyz)
                is = is_idx(kxkyz_lo, ikxkyz)

                g0 = spread((mu(:) * mu(:) * maxwell_mu(ia, iz, :, is) * aj1v(:, ikxkyz) * aj1v(:, ikxkyz)), 1, nvpa) * spread(maxwell_vpa(:, is), 2, nmu) * maxwell_fac(is) &
                * ( 1.0 - 0.5 * neo_mu_fac_global(iz, :, :, is, 1) / bmag(ia, iz) )

                wgt = 8.0 * beta * spec(is)%dens_psi0 * spec(is)%temp
                call integrate_vmu(g0, iz, tmp)
                denominator_fields_neo_33(iky, ikx, iz) = denominator_fields_neo_33(iky, ikx, iz) + tmp * wgt
            end do

            call sum_allreduce(denominator_fields_neo_33)

            denominator_fields_neo_33 = 1.0 + denominator_fields_neo_33

            ! Deallocate temporary array.
            deallocate (g0)
        end if 
    end subroutine init_neo_electromagnetic_fields  


! =================================================================================================================================================================================== !
! --------------------------------------------------- Allocate the arrays needed for higher order electrostatic field solve calculations. ------------------------------------------- !
! =================================================================================================================================================================================== !

    subroutine allocate_neo_electrostatic_fields
        use grids_z, only: nzgrid
        use grids_kxky, only: naky, nakx        

        ! Arrays. 
        use arrays, only: denominator_fields_neo

        implicit none

        if (.not. allocated(denominator_fields_neo)) then; allocate (denominator_fields_neo(naky, nakx, -nzgrid:nzgrid)); denominator_fields_neo = 0.; end if 
    end subroutine allocate_neo_electrostatic_fields

! =================================================================================================================================================================================== !
! ------------------------------------------------- Allocate the arrays needed for higher order electromagnetic field solve calculations. ------------------------------------------- !
! =================================================================================================================================================================================== !

    subroutine allocate_neo_electromagnetic_fields
        ! Grids.
        use grids_z, only: nzgrid
        use grids_kxky, only: naky, nakx
      
        ! Parameters. 
        use parameters_physics, only: include_apar, include_bpar

        ! Arrays. 
        use arrays_fields, only: apar, apar_old
        use arrays, only: denominator_fields_neo, denominator_fields_neo_12, denominator_fields_neo_13
        use arrays, only: denominator_fields_neo_21, denominator_fields_neo_22_g, denominator_fields_neo_22_gbar, denominator_fields_neo_23
        use arrays, only: denominator_fields_neo_31, denominator_fields_neo_32, denominator_fields_neo_33
  
        implicit none

        if (include_apar) then
            if (.not. allocated(denominator_fields_neo_12)) then; allocate (denominator_fields_neo_12(naky, nakx, -nzgrid:nzgrid)); denominator_fields_neo_12 = 0. ; end if
            if (.not. allocated(denominator_fields_neo_21)) then; allocate (denominator_fields_neo_21(naky, nakx, -nzgrid:nzgrid)); denominator_fields_neo_21 = 0. ; end if 
            if (.not. allocated(denominator_fields_neo_22_g)) then; allocate (denominator_fields_neo_22_g(naky, nakx, -nzgrid:nzgrid)); denominator_fields_neo_22_g = 0. ; end if
            if (.not. allocated(denominator_fields_neo_22_gbar)) then; allocate (denominator_fields_neo_22_gbar(naky, nakx, -nzgrid:nzgrid)); denominator_fields_neo_22_gbar = 0. ; end if 
        end if        
   
        if (include_bpar) then
            if (.not. allocated(denominator_fields_neo_13)) then; allocate (denominator_fields_neo_13(naky, nakx, -nzgrid:nzgrid)); denominator_fields_neo_13 = 0. ; end if
            if (.not. allocated(denominator_fields_neo_31)) then; allocate (denominator_fields_neo_31(naky, nakx, -nzgrid:nzgrid)); denominator_fields_neo_31 = 0. ; end if
            if (.not. allocated(denominator_fields_neo_23)) then; allocate (denominator_fields_neo_23(naky, nakx, -nzgrid:nzgrid)); denominator_fields_neo_23 = 0. ; end if
            if (.not. allocated(denominator_fields_neo_32)) then; allocate (denominator_fields_neo_32(naky, nakx, -nzgrid:nzgrid)); denominator_fields_neo_32 = 0. ; end if
            if (.not. allocated(denominator_fields_neo_33)) then; allocate (denominator_fields_neo_33(naky, nakx, -nzgrid:nzgrid)); denominator_fields_neo_33 = 0. ; end if
        end if
    end subroutine allocate_neo_electromagnetic_fields


! =================================================================================================================================================================================== !
! ------------------------------------------------- Deallocate the arrays needed for higher order electrostatic field solve calculations. ------------------------------------------- !
! =================================================================================================================================================================================== !

    subroutine finish_neo_electrostatic_fields
        ! Arrays.
        use arrays, only: denominator_fields_neo
        
        implicit none

        if (allocated(denominator_fields_neo)) deallocate(denominator_fields_neo)
    end subroutine finish_neo_electrostatic_fields


! =================================================================================================================================================================================== !
! ----------------------------------------------- Deallocate the arrays needed for higher order electromagnetic field solve calculations. ------------------------------------------- !
! =================================================================================================================================================================================== !

    subroutine finish_neo_electromagnetic_fields
        ! Arrays.
        use arrays, only: denominator_fields_neo, denominator_fields_neo_12, denominator_fields_neo_13
        use arrays, only: denominator_fields_neo_21, denominator_fields_neo_22_g, denominator_fields_neo_22_gbar, denominator_fields_neo_23
        use arrays, only: denominator_fields_neo_31, denominator_fields_neo_32, denominator_fields_neo_33

        implicit none

        if (allocated(denominator_fields_neo)) deallocate(denominator_fields_neo)
        if (allocated(denominator_fields_neo_12)) deallocate(denominator_fields_neo_12)
        if (allocated(denominator_fields_neo_13)) deallocate(denominator_fields_neo_13)
        if (allocated(denominator_fields_neo_21)) deallocate(denominator_fields_neo_21)
        if (allocated(denominator_fields_neo_22_g)) deallocate(denominator_fields_neo_22_g)
        if (allocated(denominator_fields_neo_22_gbar)) deallocate(denominator_fields_neo_22_gbar)
        if (allocated(denominator_fields_neo_23)) deallocate(denominator_fields_neo_23)
        if (allocated(denominator_fields_neo_31)) deallocate(denominator_fields_neo_31)
        if (allocated(denominator_fields_neo_32)) deallocate(denominator_fields_neo_32)
        if (allocated(denominator_fields_neo_33)) deallocate(denominator_fields_neo_33)
    end subroutine finish_neo_electromagnetic_fields

! =================================================================================================================================================================================== !

end module field_equations_fluxtube_neoclassical
