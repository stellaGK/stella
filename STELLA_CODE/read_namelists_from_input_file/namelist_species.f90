!###############################################################################
!################## READ STELLA NAMELISTS FOR KINETIC SPECIES ##################
!###############################################################################
! 
!                                 INPUT FILE
! 
! This module will read the namelists associated with the kinetic species:
! 
!   species_options
!     nspec = 2.0
!     species_option = 'stella'
!     read_profile_variation = .false.
!     write_profile_variation = .false.
!   
!   species_parameters_1
!     z = 1.0
!     mass = 1.0
!     dens = 1.0
!     temp = 1.0
!     tprim = -999.9
!     fprim = -999.9
!     d2ndr2 = 0.0
!     d2tdr2 = 0.0
!     bess_fac = 1.0
!     type = 'default'
!   
!   species_parameters_2
!     z = 1.0
!     mass = 1.0
!     dens = 1.0
!     temp = 1.0
!     tprim = -999.9
!     fprim = -999.9
!     d2ndr2 = 0.0
!     d2tdr2 = 0.0
!     bess_fac = 1.0
!     type = 'default'
!   
!   euterpe_parameters
!     nradii = 1000.0
!     data_file = 'euterpe.dat'
! 
!   adiabatic_electron_response
!     adiabatic_option = 'field-line-average-term'
!     tite = 1.0
!     nine = 1.0
!   
!   adiabatic_ion_response
!     adiabatic_option = 'no-field-line-average-term'
!     tite = 1.0
!     nine = 1.0
! 
!                             TEXT OPTIONS
!   
! Text options for <initialise_distribution_option>:
!    - Ions: {default, ion}
!    - Electrons: {electron, e}
! 
! Text options for <adiabatic_option>:
!    - Adiabatic electron response: {default, field-line-average-term, iphi00=2}
!    - Adiabatic ion response: {no-field-line-average-term, iphi00=0, iphi00=1}
! 
!                          ADIABATIC SPECIES
! 
! Adiabatic options: This is used when no kinetic electrons are present. The 
! non-kineticspecies (usually electrons) is set to have an adiabatic response.
! This can be either the classic adiabatic option, or the modified
! adiabatic option (i.e. modified Boltzmann electrons).
! 
!                              MODUlE
! 
! For each namelists two (or three) routines exist:
!    - set_default_parameters_<namelist>
!    - read_namelist_<namelist>
!    - check_inputs_<namelist>
! 
! First the default input parameters are set, then the default options are
! overwritten with those specified in the input file. Optionally, it is
! checked whether any input variables are clashing with each other.
!###############################################################################
module namelist_species
   
   use stella_common_types, only: spec_type

   implicit none

   ! Make reading routines accesible to other modules
   public :: read_namelist_species_options
   public :: read_namelist_species_stella
   public :: read_namelist_euterpe_parameters
   public :: read_namelist_adiabatic_electron_response
   public :: read_namelist_adiabatic_ion_response
   
   ! Parameters need to be public (species_option)
   public :: species_option_stella, species_option_inputprofs, species_option_euterpe, species_option_multibox
   
   ! Parameters need to be public (type)
   public :: ion_species, electron_species
   
   ! Parameters need to be public (adiabatic_option)
   public :: adiabatic_option_periodic, adiabatic_option_fieldlineavg

   private

   ! Create parameters for <species_option>
   integer, parameter :: species_option_stella = 1
   integer, parameter :: species_option_inputprofs = 2
   integer, parameter :: species_option_euterpe = 3
   integer, parameter :: species_option_multibox = 4

   ! Create parameters for <type>
   integer, parameter :: ion_species = 1
   integer, parameter :: electron_species = 2

   ! Create parameters for <adiabatic_option>
   integer, parameter :: adiabatic_option_periodic = 1
   integer, parameter :: adiabatic_option_fieldlineavg = 2
   
   ! These variables are used in every single subroutine, so make them global
   integer :: in_file
   logical :: dexist

contains

   !****************************************************************************
   !                              SPECIES OPTIONS                              !
   !****************************************************************************
   subroutine read_namelist_species_options(nspec, species_option_switch, &
      read_profile_variation, write_profile_variation)

      use mp, only: proc0, mp_abort
      use parameters_physics, only: initialised_parameters_physics

      implicit none

      ! Variables that are read from the input file
      integer, intent (out) :: nspec
      logical, intent (out) :: read_profile_variation, write_profile_variation
      integer, intent (out) :: species_option_switch
      
      ! Local variable to set <species_option_switch>
      character(20) :: species_option
      
      !-------------------------------------------------------------------------

      ! Some species flags are based on <radial_variation>
      ! Therefore, we need to read the physics parameters first
      if (.not. initialised_parameters_physics) then
         call mp_abort('Initialise physics parameters before reading species namelists. Aborting.')
      end if
      
      ! Set the default input parameters and read the input file
      if (.not. proc0) return
      call set_default_species_options
      call read_input_file_species_options
      call check_inputs_species_options
   
   contains
   
      !------------------------ Default input parameters -----------------------
      subroutine set_default_species_options

         implicit none

         nspec = 2
         species_option = 'stella'
         read_profile_variation = .false.
         write_profile_variation = .false.

      end subroutine set_default_species_options

      !------------------------ Read input file species options ----------------
      subroutine read_input_file_species_options

         use file_utils, only: input_unit_exist, error_unit
         use text_options, only: text_option, get_option_value

         implicit none

         ! Variables needed to read the input file
         integer :: ierr

         ! Link text options for <species_option> to an integer value
         type(text_option), dimension(4), parameter :: species_opts = (/ &
              text_option('default', species_option_stella), &
              text_option('stella', species_option_stella), &
              text_option('input.profiles', species_option_inputprofs), &
              text_option('euterpe', species_option_euterpe)/)

         ! Variables in the <species_options> namelist
         namelist /species_options/ nspec, read_profile_variation, &
            write_profile_variation, species_option
            
         !----------------------------------------------------------------------

         ! Overwrite the default input parameters by those specified in the input file
         in_file = input_unit_exist('species_options', dexist)
         if (dexist) read (unit=in_file, nml=species_options)
         
         ! Read the text option in <species_option> and store it in <species_option_switch>
         ierr = error_unit()
         call get_option_value(species_option, species_opts, species_option_switch, &
               ierr, 'species_option in namelist_species_options')

      end subroutine read_input_file_species_options

      !------------------------ Check inputs species options -------------------
      subroutine check_inputs_species_options

         use file_utils, only: runtype_option_switch, runtype_multibox
         use file_utils, only: error_unit
         use mp, only: mp_abort, job
         use parameters_physics, only: radial_variation

         implicit none

         ! Make sure we have at least one kinetic species
         if (nspec < 1) call mp_abort('Invalid nspec in species_options. Aborting.')
         
         ! Will need to readjust the species parameters in the left/right boxes
         if (runtype_option_switch == runtype_multibox .and. (job /= 1) .and. radial_variation) then
               species_option_switch = species_option_multibox
         end if

      end subroutine check_inputs_species_options

   end subroutine read_namelist_species_options
   
   !****************************************************************************
   !                        ADIABATIC ELECTRON RESPONSE                        !
   !****************************************************************************
   subroutine read_namelist_adiabatic_electron_response(adiabatic_option_switch, tite, nine)

      use mp, only: proc0

      implicit none

      ! Variables that are read from the input file
      integer, intent(out) :: adiabatic_option_switch
      real, intent(out) :: tite, nine

      ! Local variable to set <adiabatic_option_switch>
      character(30):: adiabatic_option
      
      !-------------------------------------------------------------------------

      if (.not. proc0) return
      call set_default_parameters_adiabatic_electron_response
      call read_input_file_adiabatic_electron_response
      call check_adiabatic_electron_response

   contains
      
      !------------------------ Default input parameters -----------------------
      subroutine set_default_parameters_adiabatic_electron_response

         implicit none

         ! If no 'electron' species are specified in the kinetic species, then
         ! adiabatic electrons will be added. They need a modified Boltzmann response.
         ! The electron density and temperature are set through Ti/Te and ni/ne.
         adiabatic_option = 'field-line-average-term'
         tite = 1.0
         nine = 1.0

      end subroutine set_default_parameters_adiabatic_electron_response

      !---------------------------- Read input file ----------------------------
      subroutine read_input_file_adiabatic_electron_response

         use file_utils, only: input_unit_exist, error_unit
         use text_options, only: text_option, get_option_value

         implicit none

         ! Variables needed to read the input file
         integer :: ierr

         ! Link text options for <adiabatic_option> to an integer value
         type(text_option), dimension(6), parameter :: adiabaticopts = &
            (/text_option('default', adiabatic_option_fieldlineavg), &
            text_option('no-field-line-average-term', adiabatic_option_periodic), &
            text_option('field-line-average-term', adiabatic_option_fieldlineavg), &
            text_option('iphi00=0', adiabatic_option_periodic), &
            text_option('iphi00=1', adiabatic_option_periodic), &
            text_option('iphi00=2', adiabatic_option_fieldlineavg)/)
            
         ! Variables in the <adiabatic_electron_response> namelist
         namelist /adiabatic_electron_response/ adiabatic_option, tite, nine

         !----------------------------------------------------------------------
         
         ! Overwrite the default input parameters by those specified in the input file
         in_file = input_unit_exist('adiabatic_electron_response', dexist)
         if (dexist) read (unit=in_file, nml=adiabatic_electron_response)

         ! Read the text option in <adiabatic_option> and store it in <adiabatic_option_switch>
         ierr = error_unit()
         call get_option_value(adiabatic_option, adiabaticopts, adiabatic_option_switch, &
            ierr, 'adiabatic_option in namelist_parameters.f90')

      end subroutine read_input_file_adiabatic_electron_response
      
      !--------------------------- Check parameters ----------------------------
      subroutine check_adiabatic_electron_response
      
         use mp, only: mp_abort
      
         implicit none
         
         if (adiabatic_option_switch /= adiabatic_option_fieldlineavg) then
            call mp_abort('If the adiabatic electron namelist is read, it means that no kinetic ion species &
               &have been specified. Therefore, adiabatic electrons are added which require a Modified &
               &Boltzmann response. Aborting.')
         end if
            
      end subroutine check_adiabatic_electron_response

   end subroutine read_namelist_adiabatic_electron_response
   
   !****************************************************************************
   !                           ADIABATIC ION RESPONSE                          !
   !****************************************************************************
   subroutine read_namelist_adiabatic_ion_response(adiabatic_option_switch, tite, nine)

      use mp, only: proc0

      implicit none

      ! Variables that are read from the input file
      integer, intent(out) :: adiabatic_option_switch
      real, intent(out) :: tite, nine

      ! Local variable to set <adiabatic_option_switch>
      character(30):: adiabatic_option
      
      !-------------------------------------------------------------------------

      if (.not. proc0) return
      call set_default_parameters_adiabatic_ion_response
      call read_input_file_adiabatic_ion_response
      call check_adiabatic_ion_response
      
   contains
      
      !------------------------ Default input parameters -----------------------
      subroutine set_default_parameters_adiabatic_ion_response

         implicit none

         ! If no 'ion' species are specified in the kinetic species, then
         ! adiabatic ions will be added. They need a normal Boltzmann response.
         ! The ion density and temperature are set through Ti/Te and ni/ne.
         adiabatic_option = 'no-field-line-average-term'
         tite = 1.0
         nine = 1.0

      end subroutine set_default_parameters_adiabatic_ion_response

      !---------------------------- Read input file ----------------------------
      subroutine read_input_file_adiabatic_ion_response

         use file_utils, only: input_unit_exist, error_unit
         use text_options, only: text_option, get_option_value

         implicit none

         ! Variables needed to read the input file
         integer :: ierr

         ! Link text options for <adiabatic_option> to an integer value
         type(text_option), dimension(6), parameter :: adiabaticopts = &
            (/text_option('default', adiabatic_option_fieldlineavg), &
            text_option('no-field-line-average-term', adiabatic_option_periodic), &
            text_option('field-line-average-term', adiabatic_option_fieldlineavg), &
            text_option('iphi00=0', adiabatic_option_periodic), &
            text_option('iphi00=1', adiabatic_option_periodic), &
            text_option('iphi00=2', adiabatic_option_fieldlineavg)/)
            
         ! Variables in the <adiabatic_electron_response> namelist
         namelist /adiabatic_ion_response/ adiabatic_option, tite, nine

         !----------------------------------------------------------------------
         
         ! Overwrite the default input parameters by those specified in the input file
         in_file = input_unit_exist('adiabatic_ion_response', dexist)
         if (dexist) read (unit=in_file, nml=adiabatic_ion_response)

         ! Read the text option in <adiabatic_option> and store it in <adiabatic_option_switch>
         ierr = error_unit()
         call get_option_value(adiabatic_option, adiabaticopts, adiabatic_option_switch, &
            ierr, 'adiabatic_option in namelist_parameters.f90')

      end subroutine read_input_file_adiabatic_ion_response
      
      !--------------------------- Check parameters ----------------------------
      subroutine check_adiabatic_ion_response
      
         use mp, only: mp_abort
      
         implicit none
         
         if (adiabatic_option_switch /= adiabatic_option_periodic) then
            call mp_abort('If the adiabatic ion namelist is read, it means that no kinetic electron species &
               &have been specified. Therefore, adiabatic ions are added which require a  &
               &periodic Boltzmann response. Aborting.')
         end if
            
      end subroutine check_adiabatic_ion_response

   end subroutine read_namelist_adiabatic_ion_response

   !****************************************************************************
   !                          SPECIES OPTION: STELLA                           !
   !****************************************************************************
   subroutine read_namelist_species_stella (nspec, spec)

      use mp, only: proc0

      implicit none

      ! Variables that are read from the input file
      integer, intent (in) :: nspec
      
      ! Note that <spec> is a derived variable of type <spec_type>
      ! which contains the following members: z, mass, dens, temp, ...
      type(spec_type), dimension(:), intent (out) :: spec
      real :: z, mass, dens, temp, tprim, fprim, d2ndr2, d2Tdr2, bess_fac
      
      ! Local variable to set <spec(is)%type>
      character(20) :: type
      
      !-------------------------------------------------------------------------
      
      if (.not. proc0) return
      call read_input_file_species_stella

   contains
      
      !------------------------ Set default species parameters ----------------
      subroutine set_default_species_parameters
      
         implicit none

         z = 1
         mass = 1.0
         dens = 1.0
         temp = 1.0
         tprim = -999.9
         fprim = -999.9
         d2ndr2 = 0.0
         d2Tdr2 = 0.0
         bess_fac = 1.0
         type = 'default'

      end subroutine set_default_species_parameters
   
      !------------------------ Read input file species ------------------------
      subroutine read_input_file_species_stella

         use file_utils, only: input_unit_exist, error_unit
         use file_utils, only: get_indexed_namelist_unit
         use text_options, only: text_option, get_option_value

         implicit none

         ! Variables needed to read the input file
         integer :: ierr
         
         ! Variables to read various species_parameters_<is> namelists
         integer :: unit, is
         
         ! Link text options for <type> to an integer value
         type(text_option), dimension(4), parameter :: typeopts = (/ &
              text_option('default', ion_species), &
              text_option('ion', ion_species), &
              text_option('electron', electron_species), &
              text_option('e', electron_species)/)

         ! Variables in the <species_parameters> namelist
         namelist /species_parameters/ z, mass, dens, temp, &
               tprim, fprim, d2ndr2, d2Tdr2, bess_fac, type

         !----------------------------------------------------------------------

         ! For each kinetic species, read it's namelist species_parameters_<is>
         do is = 1, nspec
         
            ! Set the default species variables
            call set_default_species_parameters
            
            ! Read the namelist species_parameters_<is>
            call get_indexed_namelist_unit(unit, 'species_parameters', is)
            read (unit=unit, nml=species_parameters)
            close (unit=unit)

            ! Save the values in the namelists as members of <spec_type>
            spec(is)%z = z
            spec(is)%mass = mass
            spec(is)%dens = dens
            spec(is)%temp = temp
            spec(is)%tprim = tprim
            spec(is)%fprim = fprim
            spec(is)%d2ndr2 = d2ndr2      ! this is (1/n_s)*d^2 n_s / drho^2
            spec(is)%d2Tdr2 = d2Tdr2      ! this is (1/T_s)*d^2 T_s / drho^2
            spec(is)%dens_psi0 = dens
            spec(is)%temp_psi0 = temp
            spec(is)%bess_fac = bess_fac

            ! Read the text option in <type> and store it in <spec(is)%type>
            ierr = error_unit()
            call get_option_value(type, typeopts, spec(is)%type, ierr, 'type in species_parameters_x')
               
         end do

      end subroutine read_input_file_species_stella

   end subroutine read_namelist_species_stella

   !****************************************************************************
   !                         SPECIES OPTIONS: EUTERPE                          !
   !****************************************************************************
   subroutine read_namelist_euterpe_parameters(nradii, data_file)

      use mp, only: proc0

      implicit none

      ! Variables that are read from the input file
      integer, intent(out) :: nradii
      character(*), intent(out) :: data_file
      
      !-------------------------------------------------------------------------
      
      if (.not. proc0) return
      call set_default_euterpe_parameters
      call read_input_file_euterpe_parameters
   
   contains
   
      !------------------------ Default input parameters -----------------------
      subroutine set_default_euterpe_parameters

         implicit none

         nradii = 1000
         data_file = 'euterpe.dat'

      end subroutine set_default_euterpe_parameters

      !------------------------ Read input file species options ----------------
      subroutine read_input_file_euterpe_parameters

         use file_utils, only: input_unit_exist, error_unit

         implicit none

         namelist /euterpe_parameters/ nradii, data_file
         in_file = input_unit_exist('euterpe_parameters', dexist)
         if (dexist) read (unit=in_file, nml=euterpe_parameters)

      end subroutine read_input_file_euterpe_parameters

   end subroutine read_namelist_euterpe_parameters

end module namelist_species
