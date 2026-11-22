!###############################################################################
!#################### READ STELLA NAMELISTS FOR DISSIPATION ####################
!###############################################################################
! 
! This module will read the namelists associated with dissipation:
! 
!   dissipation_and_collisions_options
!     include_collisions = .false.
!     collisions_implicit = .true.
!     hyper_dissipation = .false.
!     collision_model = 'dougherty'
!     ecoll_zeff = .false.
!     zeff = 1.0
!     vnew_ref = -1.0
!   
!   collisions_dougherty
!     momentum_conservation = .true.
!     energy_conservation = .true.
!     vpa_operator = .true.
!     mu_operator = .true.
!   
!   collisions_fokker_planck
!     include_field_particle_part = collisions_implicit
!     include_flr_terms = .true.
!     include_mixed_derivative_terms = .true.
!     exact_conservation = .false.
!   
!   collisions_fokker_planck_developers
!     use_explicit_algorithm = .false.
!     include_test_particle_part = .true.
!     field_particle_lmax = 1
!     field_particle_jmax = 1
!     include_like_species_collisions = .true.
!     include_unlike_species_collisions = .true.
!     scale_nu_ii = 1.0
!     scale_nu_ie = 1.0
!     scale_nu_ee = 1.0
!     scale_nu_ei = 1.0
!     include_pitch_angle_scattering = .true.
!     scale_pitch_angle_scattering_ei = 1.0
!     scale_energy_diffusion_ei = 1.0
!     include_vpa_derivatives = .true.
!     include_mu_derivatives = .true.
!     include_j1_l1_term = .true.
!     include_j1_l2_term = .true.
!     include_j0_l2_term = .true.
!     mass_ratio_approximation_ei = .false.
!     conservative_test_particle_density = .false.
!     exact_field_particle_density = .false.
!     field_particle_density_from_test_particle = .false.
!     spitzer_problem = .false.
!     spitzer_scale_epar = 1.0
!     spitzer_scale_grad_t = 0.0
!     scale_flr_terms = 1.0
!     scale_mixed_derivative_terms = 1.0
!     exact_conservation_algorithm = 'automatic'
!   
!   hyper_dissipation
!     d_hyper = 0.05
!     d_zed = 0.05
!     d_vpa = 0.05
!     hyp_zed = .false.
!     hyp_vpa = .false.
!     use_physical_ksqr = .true.
!     scale_to_outboard = .false.
! 
! Text options for <dissipation:collision_model> are:
!   {dougherty, fokker-planck}
! 
! For each namelists two (or three) routines exist:
!    - set_default_parameters_<namelist>
!    - read_namelist_<namelist>
!    - check_inputs_<namelist>
! 
! First the default input parameters are set, then the default options are
! overwritten with those specified in the input file. Optionally, it is
! checked whether any input variables are clashing with each other.
! 
!###############################################################################
module namelist_dissipation

   implicit none

   ! Make reading routines accesible to other modules
   public :: read_namelist_dissipation_and_collisions_options
   public :: read_namelist_collisions_dougherty
   public :: read_namelist_collisions_fokker_planck
   public :: read_namelist_hyper_dissipation
   
   private
   
   ! These variables are used in every single subroutine, so make them global
   integer :: in_file
   logical :: dexist

contains

   !****************************************************************************
   !                                DISSIPATION                                !
   !****************************************************************************
   subroutine read_namelist_dissipation_and_collisions_options(include_collisions, &
      collisions_implicit, collision_model, hyper_dissipation, ecoll_zeff, zeff, vnew_ref)

      use mp, only: proc0
      
      implicit none

      ! Variables that are read from the input file
      logical, intent(out) :: include_collisions, collisions_implicit, hyper_dissipation
      character(30), intent(out) :: collision_model
      logical, intent (out) :: ecoll_zeff
      real, intent(out) :: zeff, vnew_ref
         
      !-------------------------------------------------------------------------

      if (.not. proc0) return
      call set_default_parameters_dissipation_and_collisions_options
      call read_input_file_dissipation_and_collisions_options
      call check_inputs_dissipation_and_collisions_options
      call write_parameters_to_input_file

   contains 
      
      !------------------------ Default input parameters -----------------------
      subroutine set_default_parameters_dissipation_and_collisions_options

         implicit none

         ! By default we do not include collisions nor dissipation
         include_collisions = .false.
         collisions_implicit = .true.
         hyper_dissipation = .false.
         ecoll_zeff = .false.
         zeff = 1.0
         vnew_ref = -1.0
         
         ! Text options: dougherty or fokker-planck
         collision_model = 'dougherty'

      end subroutine set_default_parameters_dissipation_and_collisions_options

      !---------------------------- Read input file ----------------------------
      subroutine read_input_file_dissipation_and_collisions_options

         use file_utils, only: input_unit_exist
         implicit none

         namelist /dissipation_and_collisions_options/ include_collisions, collisions_implicit, &
            collision_model, hyper_dissipation, ecoll_zeff, zeff, vnew_ref
         in_file = input_unit_exist('dissipation_and_collisions_options', dexist)
         if (dexist) read (unit=in_file, nml=dissipation_and_collisions_options)

      end subroutine read_input_file_dissipation_and_collisions_options

      !------------------------- Check input parameters ------------------------
      subroutine check_inputs_dissipation_and_collisions_options

         implicit none

         if (.not. include_collisions) collisions_implicit = .false.

      end subroutine check_inputs_dissipation_and_collisions_options
      
      !------------------------- Write input parameters ------------------------
      subroutine write_parameters_to_input_file

         use file_units, only: unit => unit_input_file_with_defaults

         implicit none

         !-------------------------------------------------------------------------

         write (unit, '(A)') '&dissipation_and_collisions_options'
         write (unit, '(A, L1)') '  hyper_dissipation = ', hyper_dissipation
         write (unit, '(A, L1)') '  include_collisions = ', include_collisions
         write (unit, '(A, A, A)') '  collision_model = "', trim(collision_model),'"'
         write (unit, '(A, L1)') '  collisions_implicit = ', collisions_implicit
         write (unit, '(A, L1)') '  ecoll_zeff = ', ecoll_zeff
         write (unit, '(A, ES12.4)') '  zeff = ', zeff
         write (unit, '(A, ES12.4)') '  vnew_ref = ', vnew_ref
         write (unit, '(A)') '/'
         write (unit, '(A)') ''

      end subroutine write_parameters_to_input_file

   end subroutine read_namelist_dissipation_and_collisions_options

   !****************************************************************************
   !                           COLLISIONS: DOUGHERTY                           !
   !****************************************************************************
   subroutine read_namelist_collisions_dougherty(momentum_conservation, &
      energy_conservation, vpa_operator, mu_operator)

      use mp, only: proc0
      
      implicit none

      ! Variables that are read from the input file
      logical, intent(out) :: momentum_conservation, energy_conservation, vpa_operator, mu_operator
         
      !-------------------------------------------------------------------------

      if (.not. proc0) return
      call set_default_parameters_collisions_dougherty
      call read_input_file_collisions_dougherty
      call write_parameters_to_input_file
      
   contains
      
      !------------------------ Default input parameters -----------------------
      subroutine set_default_parameters_collisions_dougherty

         implicit none

         momentum_conservation = .true.       ! momentum conservation for Dougherty operator
         energy_conservation = .true.         ! energy conservation for Dougherty operator
         vpa_operator = .true.                ! include vpa components in Dougherty operator
         mu_operator = .true.                 ! include mu components in Dougherty operator

      end subroutine set_default_parameters_collisions_dougherty

      !---------------------------- Read input file ----------------------------
      subroutine read_input_file_collisions_dougherty

         use file_utils, only: input_unit_exist

         implicit none

         namelist /collisions_dougherty/ momentum_conservation, energy_conservation, vpa_operator, mu_operator
         in_file = input_unit_exist('collisions_dougherty', dexist)
         if (dexist) read (unit=in_file, nml=collisions_dougherty)

      end subroutine read_input_file_collisions_dougherty
      
      !------------------------- Write input parameters ------------------------
      subroutine write_parameters_to_input_file

         use file_units, only: unit => unit_input_file_with_defaults

         implicit none

         !-------------------------------------------------------------------------

         write (unit, '(A)') '&collisions_dougherty'
         write (unit, '(A, L1)') '  momentum_conservation = ', momentum_conservation
         write (unit, '(A, L1)') '  energy_conservation = ', energy_conservation
         write (unit, '(A, L1)') '  vpa_operator = ', vpa_operator
         write (unit, '(A, L1)') '  mu_operator = ', mu_operator
         write (unit, '(A)') '/'
         write (unit, '(A)') ''
         
      end subroutine write_parameters_to_input_file

   end subroutine read_namelist_collisions_dougherty

   !****************************************************************************
   !                          COLLISIONS: FOKKER-PLANCK                        !
   !****************************************************************************
   ! The options of the Fokker-Planck operator are read from two namelists:
   !   <collisions_fokker_planck> holds the choices a user makes: whether to include the field-particle part,
   !      the finite-Larmor-radius terms and the mixed-derivative terms, and whether to conserve momentum and
   !      energy to machine precision;
   !   <collisions_fokker_planck_developers> holds the options that change the operator for tests and
   !      benchmarks of the operator itself (switching off parts or collision channels, scale factors, the
   !      explicit algorithm, the Spitzer problem). stella prints a warning whenever one of them is changed.
   !
   ! Before Nov 2026 all options were in <collisions_fokker_planck>, under other names (e.g. cfac, nuxfac,
   ! iiknob, no_j1l1). These old names are still accepted in <collisions_fokker_planck>, with a note that
   ! tells which new option to use; setting both the old and the new name of an option stops stella.
   ! Fortran namelists cannot tell an option that was left out from an option that was set to its default,
   ! so an old name counts as set when its value differs from the default of the new option.
   !
   ! The values are returned under the internal names of the module collisions_fokkerplanck, which are the old
   ! names; cfac, nuxfac and deflknob combine the include_* switch with the scale_* factor of the same terms.
   !****************************************************************************
   subroutine read_namelist_collisions_fokker_planck(collisions_implicit, testpart, fieldpart, lmax, jmax, &
      interspec, intraspec, iiknob, ieknob, eeknob, eiknob, eiediffknob, eideflknob, deflknob, eimassr_approx, &
      spitzer_problem, density_conservation, density_conservation_field, density_conservation_tp, &
      exact_conservation, exact_conservation_algorithm, vpa_operator, mu_operator, &
      cfac, nuxfac, i1fac, i2fac, no_j1l1, no_j1l2, no_j0l2, explicit_algorithm)

      use mp, only: proc0, mp_abort

      implicit none

      ! Whether <dissipation_and_collisions_options> asks for implicit collisions
      logical, intent (in) :: collisions_implicit

      ! Variables that are read from the input file, under their internal (old) names
      logical, intent (out) :: fieldpart, testpart
      integer, intent (out) :: jmax, lmax
      logical, intent (out) :: interspec, intraspec
      logical, intent (out) :: eimassr_approx, spitzer_problem
      logical, intent (out) :: density_conservation, density_conservation_field, density_conservation_tp
      logical, intent (out) :: exact_conservation
      character(len=30), intent (out) :: exact_conservation_algorithm
      logical, intent (out) :: vpa_operator, mu_operator
      logical, intent (out) :: no_j1l1, no_j1l2, no_j0l2
      real, intent (out) :: iiknob, ieknob, eeknob, eiknob, eiediffknob, eideflknob, deflknob
      real, intent (out) :: cfac, nuxfac, i1fac, i2fac
      logical, intent (out) :: explicit_algorithm

      ! Options of <collisions_fokker_planck>
      logical :: include_field_particle_part, include_flr_terms, include_mixed_derivative_terms

      ! Options of <collisions_fokker_planck_developers>
      logical :: use_explicit_algorithm
      logical :: include_test_particle_part
      integer :: field_particle_lmax, field_particle_jmax
      logical :: include_like_species_collisions, include_unlike_species_collisions
      real :: scale_nu_ii, scale_nu_ie, scale_nu_ee, scale_nu_ei
      logical :: include_pitch_angle_scattering
      real :: scale_pitch_angle_scattering_ei, scale_energy_diffusion_ei
      logical :: include_vpa_derivatives, include_mu_derivatives
      logical :: include_j1_l1_term, include_j1_l2_term, include_j0_l2_term
      logical :: mass_ratio_approximation_ei
      logical :: conservative_test_particle_density, exact_field_particle_density
      logical :: field_particle_density_from_test_particle
      ! spitzer_problem of <collisions_fokker_planck_developers> (the old spitzer_problem has the same name)
      logical :: spitzer_problem_developers
      real :: spitzer_scale_epar, spitzer_scale_grad_t
      real :: scale_flr_terms, scale_mixed_derivative_terms

      ! Removed or renamed options that are still accepted in <collisions_fokker_planck>
      integer :: nvel_local
      logical :: advfield_coll, exact_conservation_tp

      ! The default of the field-particle part: on with the implicit algorithm (the only one that has it)
      logical :: default_field_particle_part
      ! Names of the developer options that were changed, for the warning
      character(len=1000) :: changed

      !-------------------------------------------------------------------------

      if (.not. proc0) return
      call set_default_parameters_collisions_fokker_planck_developers
      call read_input_file_collisions_fokker_planck_developers
      explicit_algorithm = use_explicit_algorithm .or. .not. collisions_implicit
      default_field_particle_part = .not. explicit_algorithm
      call set_default_parameters_collisions_fokker_planck
      call read_input_file_collisions_fokker_planck
      call combine_old_and_new_options
      call check_inputs_collisions_fokker_planck
      call write_parameters_to_input_file

   contains

      !------------------------ Default input parameters -----------------------
      subroutine set_default_parameters_collisions_fokker_planck

         implicit none

         ! The field-particle part restores the conservation of momentum and energy that the test-particle part
         ! alone violates; without it, zonal flows decay too fast (about 3 times faster at nu_ii = 0.01)
         include_field_particle_part = default_field_particle_part
         ! The finite-Larmor-radius terms: the classical diffusion of the test-particle part and the J_1 terms of
         ! the field-particle part; without them (.false.) the operator is the drift-kinetic one
         include_flr_terms = .true.
         ! The mixed-derivative (nu_x) terms of the test-particle part; exact_conservation needs them switched off
         include_mixed_derivative_terms = .true.
         ! Conserve momentum and energy to machine precision (needs include_mixed_derivative_terms = .false. and
         ! conservative_wgts_vpa = .true.); the algorithm is chosen from the mu grid, see collisions_fokkerplanck
         exact_conservation = .false.

         ! The old names, set to the values that mean "not set" (the defaults of the new options)
         fieldpart = default_field_particle_part
         testpart = .true.
         lmax = 1
         jmax = 1
         nvel_local = -1
         interspec = .true.
         intraspec = .true.
         iiknob = 1.
         ieknob = 1.
         eeknob = 1.
         eiknob = 1.
         eiediffknob = 1.
         eideflknob = 1.
         deflknob = 1.
         eimassr_approx = .false.
         advfield_coll = .true.
         spitzer_problem = .false.
         density_conservation = .false.
         density_conservation_field = .false.
         density_conservation_tp = .false.
         exact_conservation_tp = .false.
         vpa_operator = .true.
         mu_operator = .true.
         cfac = 1.
         nuxfac = 1.
         i1fac = 1.
         i2fac = 0.
         no_j1l1 = .false.
         no_j1l2 = .false.
         no_j0l2 = .false.

      end subroutine set_default_parameters_collisions_fokker_planck

      subroutine set_default_parameters_collisions_fokker_planck_developers

         implicit none

         ! Use the explicit algorithm, which has only the test-particle part and only like-species collisions;
         ! it is kept to compare parts of the operator with the implicit algorithm
         use_explicit_algorithm = .false.
         ! The test-particle part; .false. leaves only the field-particle part, to test it alone
         include_test_particle_part = .true.
         ! The maximum l (Legendre) and j (Laguerre) of the expansion of the field-particle part
         field_particle_lmax = 1
         field_particle_jmax = 1
         ! Collisions between particles of the same species, and between different species
         include_like_species_collisions = .true.
         include_unlike_species_collisions = .true.
         ! Factors on the ion-ion, ion-electron, electron-electron and electron-ion collision frequencies
         scale_nu_ii = 1.
         scale_nu_ie = 1.
         scale_nu_ee = 1.
         scale_nu_ei = 1.
         ! Pitch-angle scattering (nu_D) of the test-particle part
         include_pitch_angle_scattering = .true.
         ! Factors on the electron-ion pitch-angle scattering and energy diffusion
         scale_pitch_angle_scattering_ei = 1.
         scale_energy_diffusion_ei = 1.
         ! The vpa and mu derivatives of the test-particle part
         include_vpa_derivatives = .true.
         include_mu_derivatives = .true.
         ! The (j, l) = (1, 1), (1, 2) and (0, 2) terms of the field-particle part; (1, 1) is needed with the
         ! finite-Larmor-radius terms to agree with GS2 and CGYRO
         include_j1_l1_term = .true.
         include_j1_l2_term = .true.
         include_j0_l2_term = .true.
         ! The mass-ratio approximation of the electron-ion pitch-angle scattering (beta)
         mass_ratio_approximation_ei = .false.
         ! Discretise the test-particle part so that it conserves density to machine precision (needs the
         ! equally spaced mu grid and conservative_wgts_vpa); exact_conservation switches it on there
         conservative_test_particle_density = .false.
         ! Make the field-particle part conserve density to machine precision (lmax and jmax < 2)
         exact_field_particle_density = .false.
         ! Add a term to the field-particle part, computed with the test-particle matrix, that conserves density,
         ! also on the Gauss-Laguerre mu grid
         field_particle_density_from_test_particle = .false.
         ! Solve the Spitzer problem, a test of the operator, with factors on E_par and grad T
         spitzer_problem_developers = .false.
         spitzer_scale_epar = 1.
         spitzer_scale_grad_t = 0.
         ! Factors on the finite-Larmor-radius and mixed-derivative terms (when they are included)
         scale_flr_terms = 1.
         scale_mixed_derivative_terms = 1.
         ! The algorithm of exact_conservation: 'automatic' chooses from the mu grid, 'block_products' is the one
         ! of the equally spaced grid, 'test_particle_matrix' the one of the Gauss-Laguerre grid (which also
         ! works on the equally spaced grid, to compare the two)
         exact_conservation_algorithm = 'automatic'

      end subroutine set_default_parameters_collisions_fokker_planck_developers

      !---------------------------- Read input file ----------------------------
      subroutine read_input_file_collisions_fokker_planck

         use file_utils, only: input_unit_exist

         implicit none

         ! Variables in the <collisions_fokker_planck> namelist: the new options, and the old names
         namelist /collisions_fokker_planck/ include_field_particle_part, include_flr_terms, &
            include_mixed_derivative_terms, exact_conservation, &
            testpart, fieldpart, lmax, jmax, nvel_local, &
            interspec, intraspec, iiknob, ieknob, eeknob, eiknob, eiediffknob, eideflknob, &
            deflknob, eimassr_approx, advfield_coll, spitzer_problem, density_conservation, &
            density_conservation_field, density_conservation_tp, exact_conservation_tp, &
            vpa_operator, mu_operator, cfac, nuxfac, i1fac, i2fac, no_j1l1, no_j1l2, no_j0l2

         !----------------------------------------------------------------------

         ! Overwrite the default input parameters by those specified in the input file
         in_file = input_unit_exist('collisions_fokker_planck', dexist)
         if (dexist) read (unit=in_file, nml=collisions_fokker_planck)

      end subroutine read_input_file_collisions_fokker_planck

      subroutine read_input_file_collisions_fokker_planck_developers

         use file_utils, only: input_unit_exist

         implicit none

         ! This spitzer_problem is the one of the developer namelist; the host's spitzer_problem is the old name
         logical :: spitzer_problem

         namelist /collisions_fokker_planck_developers/ use_explicit_algorithm, include_test_particle_part, &
            field_particle_lmax, field_particle_jmax, include_like_species_collisions, &
            include_unlike_species_collisions, scale_nu_ii, scale_nu_ie, scale_nu_ee, scale_nu_ei, &
            include_pitch_angle_scattering, scale_pitch_angle_scattering_ei, scale_energy_diffusion_ei, &
            include_vpa_derivatives, include_mu_derivatives, include_j1_l1_term, include_j1_l2_term, &
            include_j0_l2_term, mass_ratio_approximation_ei, conservative_test_particle_density, &
            exact_field_particle_density, field_particle_density_from_test_particle, spitzer_problem, &
            spitzer_scale_epar, spitzer_scale_grad_t, scale_flr_terms, scale_mixed_derivative_terms, &
            exact_conservation_algorithm

         !----------------------------------------------------------------------

         spitzer_problem = spitzer_problem_developers
         in_file = input_unit_exist('collisions_fokker_planck_developers', dexist)
         if (dexist) read (unit=in_file, nml=collisions_fokker_planck_developers)
         spitzer_problem_developers = spitzer_problem

      end subroutine read_input_file_collisions_fokker_planck_developers

      !--------------------- Combine the old and new names ---------------------
      ! On return, the internal (old) names hold the values to use
      subroutine combine_old_and_new_options

         implicit none

         character(len=*), parameter :: user = 'collisions_fokker_planck'
         character(len=*), parameter :: developers = 'collisions_fokker_planck_developers'

         !----------------------------------------------------------------------

         ! Options of <collisions_fokker_planck>
         call combine_logical('fieldpart', fieldpart, 'include_field_particle_part', include_field_particle_part, &
            default_field_particle_part, user)
         call combine_switch_and_scale('cfac', cfac, 'include_flr_terms', include_flr_terms, &
            'scale_flr_terms', scale_flr_terms)
         call combine_switch_and_scale('nuxfac', nuxfac, 'include_mixed_derivative_terms', &
            include_mixed_derivative_terms, 'scale_mixed_derivative_terms', scale_mixed_derivative_terms)
         if (exact_conservation_tp) then
            if (exact_conservation) call stop_both_set('exact_conservation_tp', 'exact_conservation')
            call note_deprecated('exact_conservation_tp', 'exact_conservation (stella chooses the algorithm '// &
               'from the mu grid)', user)
            exact_conservation = .true.
         end if

         ! Options of <collisions_fokker_planck_developers>
         call combine_logical('testpart', testpart, 'include_test_particle_part', include_test_particle_part, &
            .true., developers)
         call combine_integer('lmax', lmax, 'field_particle_lmax', field_particle_lmax, 1, developers)
         call combine_integer('jmax', jmax, 'field_particle_jmax', field_particle_jmax, 1, developers)
         call combine_logical('intraspec', intraspec, 'include_like_species_collisions', &
            include_like_species_collisions, .true., developers)
         call combine_logical('interspec', interspec, 'include_unlike_species_collisions', &
            include_unlike_species_collisions, .true., developers)
         call combine_real('iiknob', iiknob, 'scale_nu_ii', scale_nu_ii, 1., developers)
         call combine_real('ieknob', ieknob, 'scale_nu_ie', scale_nu_ie, 1., developers)
         call combine_real('eeknob', eeknob, 'scale_nu_ee', scale_nu_ee, 1., developers)
         call combine_real('eiknob', eiknob, 'scale_nu_ei', scale_nu_ei, 1., developers)
         call combine_real('eideflknob', eideflknob, 'scale_pitch_angle_scattering_ei', &
            scale_pitch_angle_scattering_ei, 1., developers)
         call combine_real('eiediffknob', eiediffknob, 'scale_energy_diffusion_ei', scale_energy_diffusion_ei, &
            1., developers)
         ! deflknob is a factor on nu_D that is only meaningful as 0 or 1
         if (deflknob /= 1.) then
            if (.not. include_pitch_angle_scattering) call stop_both_set('deflknob', 'include_pitch_angle_scattering')
            if (deflknob /= 0.) call mp_abort('deflknob must be 0 or 1 in &collisions_fokker_planck. Aborting.')
            call note_deprecated('deflknob = 0', 'include_pitch_angle_scattering = .false.', developers)
         else if (.not. include_pitch_angle_scattering) then
            deflknob = 0.
         end if
         call combine_logical('vpa_operator', vpa_operator, 'include_vpa_derivatives', include_vpa_derivatives, &
            .true., developers)
         call combine_logical('mu_operator', mu_operator, 'include_mu_derivatives', include_mu_derivatives, &
            .true., developers)
         ! no_j1l1, no_j1l2 and no_j0l2 switch the terms off, the new options switch them on
         call combine_logical('no_j1l1', no_j1l1, 'include_j1_l1_term (the opposite of no_j1l1)', &
            .not. include_j1_l1_term, .false., developers)
         call combine_logical('no_j1l2', no_j1l2, 'include_j1_l2_term (the opposite of no_j1l2)', &
            .not. include_j1_l2_term, .false., developers)
         call combine_logical('no_j0l2', no_j0l2, 'include_j0_l2_term (the opposite of no_j0l2)', &
            .not. include_j0_l2_term, .false., developers)
         call combine_logical('eimassr_approx', eimassr_approx, 'mass_ratio_approximation_ei', &
            mass_ratio_approximation_ei, .false., developers)
         call combine_logical('density_conservation', density_conservation, 'conservative_test_particle_density', &
            conservative_test_particle_density, .false., developers)
         call combine_logical('density_conservation_field', density_conservation_field, &
            'exact_field_particle_density', exact_field_particle_density, .false., developers)
         call combine_logical('density_conservation_tp', density_conservation_tp, &
            'field_particle_density_from_test_particle', field_particle_density_from_test_particle, .false., developers)
         ! spitzer_problem keeps its name, but moves to the developer namelist
         call combine_logical('spitzer_problem', spitzer_problem, 'spitzer_problem', spitzer_problem_developers, &
            .false., developers)
         call combine_real('i1fac', i1fac, 'spitzer_scale_epar', spitzer_scale_epar, 1., developers)
         call combine_real('i2fac', i2fac, 'spitzer_scale_grad_t', spitzer_scale_grad_t, 0., developers)

         ! Removed options
         if (nvel_local /= -1) then
            write (*, '(A)') 'Note: nvel_local in &collisions_fokker_planck is no longer used and is ignored '// &
               '(it set up a speed grid that was never used); remove it from the input file.'
         end if
         if (.not. advfield_coll) then
            call mp_abort('advfield_coll = .false. (the field-particle part without the electrostatic potential '// &
               'terms) has been removed from the Fokker-Planck operator. Aborting.')
         end if

      end subroutine combine_old_and_new_options

      !------------------------------ Check inputs -----------------------------
      subroutine check_inputs_collisions_fokker_planck

         implicit none

         character(len=20) :: value

         !----------------------------------------------------------------------

         if (.not. any(exact_conservation_algorithm == [character(len=30) :: 'automatic', 'block_products', &
            'test_particle_matrix'])) then
            call mp_abort('exact_conservation_algorithm in &collisions_fokker_planck_developers must be "automatic", '// &
               '"block_products" or "test_particle_matrix". Aborting.')
         end if

         ! The explicit algorithm is a developer option, chosen in <collisions_fokker_planck_developers>
         if (.not. collisions_implicit .and. .not. use_explicit_algorithm) then
            call mp_abort('The explicit Fokker-Planck algorithm (collisions_implicit = .false.) is a developer '// &
               'option: it has only the test-particle part and only like-species collisions. Use the implicit '// &
               'algorithm, or set use_explicit_algorithm = .true. in &collisions_fokker_planck_developers. Aborting.')
         end if
         if (explicit_algorithm) then
            write (*, '(A)') 'Warning: the explicit Fokker-Planck algorithm has only the test-particle part and '// &
               'only like-species collisions, so it does not conserve momentum and energy.'
         end if

         ! The field-particle part switched off with the implicit algorithm
         if (.not. explicit_algorithm .and. .not. fieldpart) then
            write (*, '(A)') 'Warning: the field-particle part of the Fokker-Planck operator is switched off '// &
               '(include_field_particle_part = .false.), so the collisions do not conserve momentum and energy; '// &
               'zonal flows then decay too fast (about 3 times faster at nu_ii = 0.01).'
         end if
         if (cfac == 0.) then
            write (*, '(A)') 'Warning: the finite-Larmor-radius terms of the Fokker-Planck operator are switched off '// &
               '(include_flr_terms = .false.): this is the drift-kinetic operator.'
         end if

         ! List the developer options that differ from their defaults
         changed = ''
         if (explicit_algorithm) call add_changed('use_explicit_algorithm', .true.)
         if (.not. testpart) call add_changed('include_test_particle_part', .false.)
         if (lmax /= 1) call add_changed_integer('field_particle_lmax', lmax)
         if (jmax /= 1) call add_changed_integer('field_particle_jmax', jmax)
         if (.not. intraspec) call add_changed('include_like_species_collisions', .false.)
         if (.not. interspec) call add_changed('include_unlike_species_collisions', .false.)
         if (iiknob /= 1.) call add_changed_real('scale_nu_ii', iiknob)
         if (ieknob /= 1.) call add_changed_real('scale_nu_ie', ieknob)
         if (eeknob /= 1.) call add_changed_real('scale_nu_ee', eeknob)
         if (eiknob /= 1.) call add_changed_real('scale_nu_ei', eiknob)
         if (deflknob /= 1.) call add_changed('include_pitch_angle_scattering', .false.)
         if (eideflknob /= 1.) call add_changed_real('scale_pitch_angle_scattering_ei', eideflknob)
         if (eiediffknob /= 1.) call add_changed_real('scale_energy_diffusion_ei', eiediffknob)
         if (.not. vpa_operator) call add_changed('include_vpa_derivatives', .false.)
         if (.not. mu_operator) call add_changed('include_mu_derivatives', .false.)
         if (no_j1l1) call add_changed('include_j1_l1_term', .false.)
         if (no_j1l2) call add_changed('include_j1_l2_term', .false.)
         if (no_j0l2) call add_changed('include_j0_l2_term', .false.)
         if (eimassr_approx) call add_changed('mass_ratio_approximation_ei', .true.)
         if (density_conservation) call add_changed('conservative_test_particle_density', .true.)
         if (density_conservation_field) call add_changed('exact_field_particle_density', .true.)
         if (density_conservation_tp) call add_changed('field_particle_density_from_test_particle', .true.)
         if (spitzer_problem) call add_changed('spitzer_problem', .true.)
         if (i1fac /= 1.) call add_changed_real('spitzer_scale_epar', i1fac)
         if (i2fac /= 0.) call add_changed_real('spitzer_scale_grad_t', i2fac)
         if (cfac /= 0. .and. cfac /= 1.) call add_changed_real('scale_flr_terms', cfac)
         if (nuxfac /= 0. .and. nuxfac /= 1.) call add_changed_real('scale_mixed_derivative_terms', nuxfac)
         if (exact_conservation_algorithm /= 'automatic') then
            value = '"'//trim(exact_conservation_algorithm)//'"'
            changed = trim(changed)//' exact_conservation_algorithm = '//trim(value)//','
         end if
         if (len_trim(changed) > 0) then
            write (*, '(A)') 'Warning: developer options of the Fokker-Planck operator changed:'// &
               changed(1:len_trim(changed) - 1)//'. These change the operator and are meant for tests.'
         end if

      end subroutine check_inputs_collisions_fokker_planck

      !------------------------- Write input parameters ------------------------
      subroutine write_parameters_to_input_file

         use file_units, only: unit => unit_input_file_with_defaults

         implicit none

         !-------------------------------------------------------------------------

         write (unit, '(A)') '&collisions_fokker_planck'
         write (unit, '(A, L1)') '  include_field_particle_part = ', fieldpart
         write (unit, '(A, L1)') '  include_flr_terms = ', cfac /= 0.
         write (unit, '(A, L1)') '  include_mixed_derivative_terms = ', nuxfac /= 0.
         write (unit, '(A, L1)') '  exact_conservation = ', exact_conservation
         write (unit, '(A)') '/'
         write (unit, '(A)') ''
         write (unit, '(A)') '&collisions_fokker_planck_developers'
         write (unit, '(A, L1)') '  use_explicit_algorithm = ', explicit_algorithm
         write (unit, '(A, L1)') '  include_test_particle_part = ', testpart
         write (unit, '(A, I0)') '  field_particle_lmax = ', lmax
         write (unit, '(A, I0)') '  field_particle_jmax = ', jmax
         write (unit, '(A, L1)') '  include_like_species_collisions = ', intraspec
         write (unit, '(A, L1)') '  include_unlike_species_collisions = ', interspec
         write (unit, '(A, ES12.4)') '  scale_nu_ii = ', iiknob
         write (unit, '(A, ES12.4)') '  scale_nu_ie = ', ieknob
         write (unit, '(A, ES12.4)') '  scale_nu_ee = ', eeknob
         write (unit, '(A, ES12.4)') '  scale_nu_ei = ', eiknob
         write (unit, '(A, L1)') '  include_pitch_angle_scattering = ', deflknob /= 0.
         write (unit, '(A, ES12.4)') '  scale_pitch_angle_scattering_ei = ', eideflknob
         write (unit, '(A, ES12.4)') '  scale_energy_diffusion_ei = ', eiediffknob
         write (unit, '(A, L1)') '  include_vpa_derivatives = ', vpa_operator
         write (unit, '(A, L1)') '  include_mu_derivatives = ', mu_operator
         write (unit, '(A, L1)') '  include_j1_l1_term = ', .not. no_j1l1
         write (unit, '(A, L1)') '  include_j1_l2_term = ', .not. no_j1l2
         write (unit, '(A, L1)') '  include_j0_l2_term = ', .not. no_j0l2
         write (unit, '(A, L1)') '  mass_ratio_approximation_ei = ', eimassr_approx
         write (unit, '(A, L1)') '  conservative_test_particle_density = ', density_conservation
         write (unit, '(A, L1)') '  exact_field_particle_density = ', density_conservation_field
         write (unit, '(A, L1)') '  field_particle_density_from_test_particle = ', density_conservation_tp
         write (unit, '(A, L1)') '  spitzer_problem = ', spitzer_problem
         write (unit, '(A, ES12.4)') '  spitzer_scale_epar = ', i1fac
         write (unit, '(A, ES12.4)') '  spitzer_scale_grad_t = ', i2fac
         write (unit, '(A, ES12.4)') '  scale_flr_terms = ', merge(cfac, 1., cfac /= 0.)
         write (unit, '(A, ES12.4)') '  scale_mixed_derivative_terms = ', merge(nuxfac, 1., nuxfac /= 0.)
         write (unit, '(A)') "  exact_conservation_algorithm = '"//trim(exact_conservation_algorithm)//"'"
         write (unit, '(A)') '/'
         write (unit, '(A)') ''

      end subroutine write_parameters_to_input_file

      !------------------------------- Helpers ---------------------------------
      ! Each combine_* routine takes the value of an old name (in old_value) and of its new name, and returns in
      ! old_value the value to use: the new one, or the old one (with a note) if only the old name was set
      subroutine combine_logical(old_name, old_value, new_name, new_value, default_value, new_namelist)

         implicit none

         character(len=*), intent(in) :: old_name, new_name, new_namelist
         logical, intent(in out) :: old_value
         logical, intent(in) :: new_value, default_value

         if (old_value .neqv. default_value) then
            if (new_value .neqv. default_value) call stop_both_set(old_name, new_name)
            call note_deprecated(old_name, new_name, new_namelist)
         else
            old_value = new_value
         end if

      end subroutine combine_logical

      subroutine combine_real(old_name, old_value, new_name, new_value, default_value, new_namelist)

         implicit none

         character(len=*), intent(in) :: old_name, new_name, new_namelist
         real, intent(in out) :: old_value
         real, intent(in) :: new_value, default_value

         if (old_value /= default_value) then
            if (new_value /= default_value) call stop_both_set(old_name, new_name)
            call note_deprecated(old_name, new_name, new_namelist)
         else
            old_value = new_value
         end if

      end subroutine combine_real

      subroutine combine_integer(old_name, old_value, new_name, new_value, default_value, new_namelist)

         implicit none

         character(len=*), intent(in) :: old_name, new_name, new_namelist
         integer, intent(in out) :: old_value
         integer, intent(in) :: new_value, default_value

         if (old_value /= default_value) then
            if (new_value /= default_value) call stop_both_set(old_name, new_name)
            call note_deprecated(old_name, new_name, new_namelist)
         else
            old_value = new_value
         end if

      end subroutine combine_integer

      ! The old factors cfac and nuxfac (default 1) became a switch in <collisions_fokker_planck> (0 or 1) and a
      ! scale factor in <collisions_fokker_planck_developers> (any other value); factor returns the combination
      subroutine combine_switch_and_scale(old_name, factor, switch_name, switch, scale_name, scale)

         implicit none

         character(len=*), intent(in) :: old_name, switch_name, scale_name
         real, intent(in out) :: factor
         logical, intent(in) :: switch
         real, intent(in) :: scale

         if (factor /= 1.) then
            if (.not. switch .or. scale /= 1.) call stop_both_set(old_name, switch_name//' or '//scale_name)
            if (factor == 0.) then
               call note_deprecated(old_name//' = 0', switch_name//' = .false.', 'collisions_fokker_planck')
            else
               call note_deprecated(old_name, scale_name, 'collisions_fokker_planck_developers')
            end if
         else if (switch) then
            factor = scale
         else
            factor = 0.
         end if

      end subroutine combine_switch_and_scale

      subroutine note_deprecated(old_name, new_name, new_namelist)

         implicit none

         character(len=*), intent(in) :: old_name, new_name, new_namelist

         write (*, '(A)') 'Note: '//old_name//' in &collisions_fokker_planck is deprecated: use '//new_name// &
            ' in &'//new_namelist//' instead.'

      end subroutine note_deprecated

      subroutine stop_both_set(old_name, new_name)

         implicit none

         character(len=*), intent(in) :: old_name, new_name

         call mp_abort('Both the deprecated option '//old_name//' and its new name '//new_name//' of the '// &
            'Fokker-Planck operator are set in the input file: set only '//new_name//'. Aborting.')

      end subroutine stop_both_set

      subroutine add_changed(name, value)

         implicit none

         character(len=*), intent(in) :: name
         logical, intent(in) :: value

         if (value) then
            changed = trim(changed)//' '//name//' = .true.,'
         else
            changed = trim(changed)//' '//name//' = .false.,'
         end if

      end subroutine add_changed

      subroutine add_changed_real(name, value)

         implicit none

         character(len=*), intent(in) :: name
         real, intent(in) :: value
         character(len=20) :: text

         write (text, '(G0.4)') value
         changed = trim(changed)//' '//name//' = '//trim(adjustl(text))//','

      end subroutine add_changed_real

      subroutine add_changed_integer(name, value)

         implicit none

         character(len=*), intent(in) :: name
         integer, intent(in) :: value
         character(len=20) :: text

         write (text, '(I0)') value
         changed = trim(changed)//' '//name//' = '//trim(text)//','

      end subroutine add_changed_integer

   end subroutine read_namelist_collisions_fokker_planck

   !****************************************************************************
   !                              HYPER DISSIPATION                            !
   !****************************************************************************
   subroutine read_namelist_hyper_dissipation(D_hyper, D_zed, D_vpa, hyp_zed, &
      hyp_vpa, use_physical_ksqr, scale_to_outboard)

      use mp, only: proc0, mp_abort
      use parameters_physics, only: initialised_parameters_physics
      
      implicit none

      ! Variables that are read from the input file
      logical, intent (out) :: use_physical_ksqr, scale_to_outboard
      real, intent (out) :: D_hyper, D_zed, D_vpa
      logical, intent (out) :: hyp_vpa, hyp_zed
         
      !-------------------------------------------------------------------------
      
      ! The <use_physical_ksqr> flag is turned off for full_flux_surface and radial_variation
      ! Therefore, we need to read the physics parameters first
      if (.not. initialised_parameters_physics) then
         call mp_abort('Initialise physics parameters before reading dissipation namelists. Aborting.')
      end if

      ! Set the default input parameters and read the input file
      if (.not. proc0) return
      call set_default_parameters_hyper_dissipation
      call read_input_file_hyper_dissipation
      call write_parameters_to_input_file
      
   contains
      
      !------------------------ Default input parameters -----------------------
      subroutine set_default_parameters_hyper_dissipation

         use parameters_physics, only: full_flux_surface, radial_variation

         implicit none

         use_physical_ksqr = .not. (full_flux_surface .or. radial_variation)  ! use kperp2, instead of akx^2 + aky^2
         scale_to_outboard = .false.                                          ! scales hyperdissipation to zed = 0
         D_hyper = 0.05
         D_zed = 0.05
         D_vpa = 0.05
         hyp_vpa = .false.
         hyp_zed = .false.

      end subroutine set_default_parameters_hyper_dissipation

      !---------------------------- Read input file ----------------------------
      subroutine read_input_file_hyper_dissipation

         use file_utils, only: input_unit_exist

         implicit none

         ! Variables in the <hyper_dissipation> namelist
         namelist /hyper_dissipation/ D_hyper, D_zed, D_vpa, hyp_zed, &
            hyp_vpa, use_physical_ksqr, scale_to_outboard
         
         !----------------------------------------------------------------------

         ! Overwrite the default input parameters by those specified in the input file
         in_file = input_unit_exist('hyper_dissipation', dexist)
         if (dexist) read (unit=in_file, nml=hyper_dissipation)

      end subroutine read_input_file_hyper_dissipation
      
      
      !------------------------- Write input parameters ------------------------
      subroutine write_parameters_to_input_file

         use file_units, only: unit => unit_input_file_with_defaults

         implicit none

         !-------------------------------------------------------------------------

         write (unit, '(A)') '&hyper_dissipation'
         write (unit, '(A, L1)') '  use_physical_ksqr = ', use_physical_ksqr
         write (unit, '(A, L1)') '  scale_to_outboard = ', scale_to_outboard
         write (unit, '(A, L1)') '  hyp_vpa = ', hyp_vpa
         write (unit, '(A, L1)') '  hyp_zed = ', hyp_zed
         write (unit, '(A, ES12.4)') '  D_hyper = ', D_hyper
         write (unit, '(A, ES12.4)') '  D_zed = ', D_zed
         write (unit, '(A, ES12.4)') '  D_vpa = ', D_vpa
         write (unit, '(A)') '/'
         write (unit, '(A)') ''

      end subroutine write_parameters_to_input_file

   end subroutine read_namelist_hyper_dissipation

end module namelist_dissipation


