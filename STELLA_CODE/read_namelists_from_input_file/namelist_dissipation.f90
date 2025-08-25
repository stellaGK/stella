!###############################################################################
!#################### READ STELLA NAMELISTS FOR DISSIPATION ####################
!###############################################################################
! 
! This module will read the namelists associated with dissipation:
! 
!   dissipation_and_collisions_options
!      include_collisions = .false.
!      collisions_implicit = .true.
!      hyper_dissipation = .false.
!      collision_model = 'dougherty'
!   
!   collisions_dougherty
!     momentum_conservation = .true.
!     energy_conservation = .true.
!     vpa_operator = .true.
!     mu_operator = .true.
!   
!   collisions_fokker_planck
!     testpart = .true.
!     fieldpart = .false.
!     lmax = 1.0
!     jmax = 1.0
!     nvel_local = 512.0
!     interspec = .true.
!     intraspec = .true.
!     iiknob = 1.0
!     ieknob = 1.0
!     eeknob = 1.0
!     eiknob = 1.0
!     eiediffknob = 1.0
!     eideflknob = 1.0
!     deflknob = 1.0
!     eimassr_approx = .false.
!     advfield_coll = .true.
!     spitzer_problem = .false.
!     density_conservation = .false.
!     density_conservation_field = .false.
!     density_conservation_tp = .false.
!     exact_conservation = .false.
!     exact_conservation_tp = .false.
!     vpa_operator = .true.
!     mu_operator = .true.
!     cfac = 1.0
!     cfac2 = 1.0
!     nuxfac = 1.0
!     i1fac = 1.0
!     i2fac = 0.0
!     no_j1l1 = .true.
!     no_j1l2 = .false.
!     no_j0l2 = .false.
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
      collisions_implicit, collision_model, hyper_dissipation)

      use mp, only: proc0
      
      implicit none

      ! Variables that are read from the input file
      logical, intent(out) :: include_collisions, collisions_implicit, hyper_dissipation
      character(30), intent(out) :: collision_model
         
      !-------------------------------------------------------------------------

      if (.not. proc0) return
      call set_default_parameters_dissipation_and_collisions_options
      call read_input_file_dissipation_and_collisions_options
      call check_inputs_dissipation_and_collisions_options

   contains 
      
      !------------------------ Default input parameters -----------------------
      subroutine set_default_parameters_dissipation_and_collisions_options

         implicit none

         ! By default we do not include collisions nor dissipation
         include_collisions = .false.
         collisions_implicit = .true.
         hyper_dissipation = .false.
         
         ! Text options: dougherty or fokker-planck
         collision_model = 'dougherty'

      end subroutine set_default_parameters_dissipation_and_collisions_options

      !---------------------------- Read input file ----------------------------
      subroutine read_input_file_dissipation_and_collisions_options

         use file_utils, only: input_unit_exist
         implicit none

         namelist /dissipation_and_collisions_options/ include_collisions, collisions_implicit, collision_model, hyper_dissipation
         in_file = input_unit_exist('dissipation_and_collisions_options', dexist)
         if (dexist) read (unit=in_file, nml=dissipation_and_collisions_options)

      end subroutine read_input_file_dissipation_and_collisions_options

      !------------------------- Check input parameters ------------------------
      subroutine check_inputs_dissipation_and_collisions_options

         implicit none

         if (.not. include_collisions) collisions_implicit = .false.

      end subroutine check_inputs_dissipation_and_collisions_options

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

   end subroutine read_namelist_collisions_dougherty

   !****************************************************************************
   !                          COLLISIONS: FOKKER-PLANCK                        !
   !****************************************************************************
   subroutine read_namelist_collisions_fokker_planck(testpart, fieldpart, lmax, jmax, nvel_local, &
      interspec, intraspec, iiknob, ieknob, eeknob, eiknob, eiediffknob, eideflknob, deflknob, eimassr_approx, &
      advfield_coll, spitzer_problem, density_conservation, density_conservation_field, density_conservation_tp, &
      exact_conservation, exact_conservation_tp, vpa_operator, mu_operator, &
      cfac, cfac2, nuxfac, i1fac, i2fac, no_j1l1, no_j1l2, no_j0l2)

      use mp, only: proc0
      
      implicit none
      
      ! Variables that are read from the input file
      logical, intent (out) :: fieldpart, testpart
      integer, intent (out) :: jmax, lmax, nvel_local
      logical, intent (out) :: interspec, intraspec
      logical, intent (out) :: eimassr_approx, advfield_coll, spitzer_problem
      logical, intent (out) :: density_conservation, density_conservation_field, density_conservation_tp
      logical, intent (out) ::exact_conservation_tp, exact_conservation
      logical, intent (out) :: vpa_operator, mu_operator
      logical, intent (out) :: no_j1l1, no_j1l2, no_j0l2
      real, intent (out) :: iiknob, ieknob, eeknob, eiknob, eiediffknob, eideflknob, deflknob
      real, intent (out) :: cfac, cfac2, nuxfac, i1fac, i2fac
         
      !-------------------------------------------------------------------------

      if (.not. proc0) return
      call set_default_parameters_collisions_fokker_planck
      call read_input_file_collisions_fokker_planck

   contains 
      
      !------------------------ Default input parameters -----------------------
      subroutine set_default_parameters_collisions_fokker_planck

         implicit none

         testpart = .true.                    ! test particle component (TPO) of fokker-planck operator, must be True
         fieldpart = .false.                  ! enable the field particle component (FPO) of the fokker-planck operator
         intraspec = .true.                   ! intra-species collisions in the Fokker-Planck operator
         interspec = .true.                   ! inter-species
         iiknob = 1.                          ! control the ion-ion coll freq in Fokker-Planck operator
         ieknob = 1.                          ! ...ion-eon coll freq
         eeknob = 1.                          ! ...eon-eon coll freq
         eiknob = 1.                          ! ...eon-ion coll freq
         eiediffknob = 1.                     ! control the eon-ion energy diffusion in Fokker-Planck operator
         eideflknob = 1.                      ! 
         deflknob = 1.                        ! control pitch angle scattering in Fokker-Planck operator, must be 1 or 0
         eimassr_approx = .false.             ! use mass ratio approximation for test particle operator, beta
         advfield_coll = .true.               ! disable electrostatic potential terms in the field particle operator, beta
         density_conservation = .false.       ! if True and equally_spaced_mu_grid=True and conservative_wgts_vpa=True, then TPO conserves density to machine precision
         density_conservation_field = .false. ! if True and jmax, lmax < 2, then FPO conserves density to machine precision
         density_conservation_tp = .false.    ! if True add term to field particle operator to ensure density conservation, also on non-uniform grids
         exact_conservation = .false.         ! if True and fieldpart=True and lmax=jmax=1 then momentum and energy conserved to machine precision - in beta &
         ! & works only if nux = 0, need to correct the discretisation of nux terms in TPO
         exact_conservation_tp = .false.      ! if True and lmax=jmax=1 then momentum and energy conserved to machine precision, by using the test particle operator &
         ! to compute field particle terms; this is slower than exact_conservation
         spitzer_problem = .false.            ! to solve the Spitzer problem for tests of the collision operator
         cfac = 1                             ! scale gyrodiffusive term in test particle component of Fokker-Planck operator
         cfac2 = 1                            ! scale gyrodiffusive terms in field particle component of Fokker-Planck operator - in beta
         nuxfac = 1                           ! scale nux (mixed derivative) terms in test particle component of Fokker-Planck operator
         jmax = 1                             ! maximum j in Hirshman-Sigmar expansion of the field particle operator
         lmax = 1                             ! maximum l in spherical harmonic expansion of the field particle operator
         i1fac = 1                            ! for Spitzer problem
         i2fac = 0                            ! for Spitzer problem
         no_j1l1 = .true.                     ! disable j1l1 term in the field particle component of Fokker-Planck operator
         no_j1l2 = .false.                    ! disable j1l2 term
         no_j0l2 = .false.                    ! disable j0l2 term
         vpa_operator = .true.                ! include vpa components in Dougherty or Fokker-Planck operator
         mu_operator = .true.                 ! include mu components in Dougherty or Fokker-Planck operator
         nvel_local = 512
         
      end subroutine set_default_parameters_collisions_fokker_planck

      !---------------------------- Read input file ----------------------------
      subroutine read_input_file_collisions_fokker_planck

         use file_utils, only: input_unit_exist

         implicit none

         ! Variables in the <collisions_fokker_planck> namelist
         namelist /collisions_fokker_planck/ testpart, fieldpart, lmax, jmax, nvel_local, &
            interspec, intraspec, iiknob, ieknob, eeknob, eiknob, eiediffknob, eideflknob, &
            deflknob, eimassr_approx, advfield_coll, spitzer_problem, density_conservation, &
            density_conservation_field, density_conservation_tp, exact_conservation, exact_conservation_tp, &
            vpa_operator, mu_operator, cfac, cfac2, nuxfac, i1fac, i2fac, no_j1l1, no_j1l2, no_j0l2
         
         !----------------------------------------------------------------------

         ! Overwrite the default input parameters by those specified in the input file
         in_file = input_unit_exist('collisions_fokker_planck', dexist)
         if (dexist) read (unit=in_file, nml=collisions_fokker_planck)

      end subroutine read_input_file_collisions_fokker_planck

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

   end subroutine read_namelist_hyper_dissipation

end module namelist_dissipation


