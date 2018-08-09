module species

  use common_types, only: spec_type

  implicit none

  public :: init_species, finish_species
!  public :: reinit_species, init_trin_species
  public :: nspec, spec
  public :: ion_species, electron_species, slowing_down_species, tracer_species
  public :: has_electron_species, has_slowing_down_species
  public :: ions, electrons, impurity

  private

  integer, parameter :: ion_species = 1
  integer, parameter :: electron_species = 2 ! for collision operator
  integer, parameter :: slowing_down_species = 3 ! slowing-down distn
  integer, parameter :: tracer_species = 4 ! for test particle diffusion studies

  integer :: species_option_switch
  integer, parameter :: species_option_stella = 1
  integer, parameter :: species_option_inputprofs = 2
  integer, parameter :: species_option_euterpe = 3

  integer :: nspec
  type (spec_type), dimension (:), allocatable :: spec

  integer :: ions, electrons, impurity
  integer :: ntspec_trin
  real, dimension (:), allocatable :: dens_trin, temp_trin, fprim_trin, tprim_trin, nu_trin

  character (20) :: species_option

  logical :: initialized = .false.

contains

  subroutine init_species

!    use mp, only: trin_flag
    use mp, only: proc0, broadcast
    use physics_parameters, only: vnew_ref, zeff
    use inputprofiles_interface, only: read_inputprof_spec
    use euterpe_interface, only: read_species_euterpe

    implicit none

    integer :: is

    if (initialized) return
    initialized = .true.

    if (proc0) call read_species_knobs
    call broadcast (nspec)
    allocate (spec(nspec))
    if (proc0) then
       select case (species_option_switch)
       case (species_option_stella)
          call read_species_stella
       case (species_option_inputprofs)
          call read_species_stella
          call read_inputprof_spec (nspec, spec)
       case (species_option_euterpe)
          call read_species_stella
          call read_species_euterpe (nspec, spec)
       end select

       do is = 1, nspec
          ! initialize nu_ss' = 0 for all s'
          spec(is)%vnew = 0.
          ! FLAG -- only contains self-collisions at the moment
          spec(is)%vnew(is) = vnew_ref*spec(is)%dens*spec(is)%z**4 &
               / (sqrt(spec(is)%mass)*spec(is)%temp**1.5)
          ! include electron-ion collisions
          if (spec(is)%type == electron_species) then
             spec(is)%vnew(is) = spec(is)%vnew(is)*(1.+zeff)
          end if
       end do

       open (1003,file='species.input',status='unknown')
       write (1003,'(7a12,a9)') '#1.z', '2.mass', '3.dens', &
            '4.temp', '5.tprim','6.fprim', '7.vnewss', '8.type'
       do is = 1, nspec
          write (1003,'(7e12.4,i9)') spec(is)%z, spec(is)%mass, &
               spec(is)%dens, spec(is)%temp, spec(is)%tprim, &
               spec(is)%fprim, spec(is)%vnew(is), spec(is)%type
       end do
       close (1003)
    end if

    call broadcast_parameters

!    if (trin_flag) call reinit_species (ntspec_trin, dens_trin, &
!         temp_trin, fprim_trin, tprim_trin, nu_trin)
  end subroutine init_species

  subroutine read_species_knobs

    use file_utils, only: error_unit, input_unit_exist
    use text_options, only: text_option, get_option_value

    implicit none

    integer :: ierr, in_file
    logical :: exist

    namelist /species_knobs/ nspec, species_option

    type (text_option), dimension (4), parameter :: specopts = (/ &
         text_option('default', species_option_stella), &
         text_option('stella', species_option_stella), &
         text_option('input.profiles', species_option_inputprofs), &
         text_option('euterpe', species_option_euterpe) /)

    nspec = 2
    species_option = 'stella'
    
    in_file = input_unit_exist("species_knobs", exist)
    if (exist) read (unit=in_file, nml=species_knobs)

    ierr = error_unit()
    call get_option_value (species_option, specopts, species_option_switch, &
         ierr, "species_option in species_knobs")

    if (nspec < 1) then
       ierr = error_unit()
       write (unit=ierr, &
            fmt="('Invalid nspec in species_knobs: ', i5)") nspec
       stop
    end if

  end subroutine read_species_knobs

  subroutine read_species_stella

    use file_utils, only: error_unit, get_indexed_namelist_unit
    use text_options, only: text_option, get_option_value

    implicit none

    real :: z, mass, dens, temp, tprim, fprim, d2ndr2, d2Tdr2
    integer :: ierr, unit, is

    namelist /species_parameters/ z, mass, dens, temp, &
         tprim, fprim, d2ndr2, d2Tdr2, type

    character(20) :: type
    type (text_option), dimension (9), parameter :: typeopts = (/ &
         text_option('default', ion_species), &
         text_option('ion', ion_species), &
         text_option('electron', electron_species), &
         text_option('e', electron_species), &
         text_option('beam', slowing_down_species), &
         text_option('fast', slowing_down_species), &
         text_option('alpha', slowing_down_species), &
         text_option('slowing-down', slowing_down_species), &
         text_option('trace', tracer_species) /)

    do is = 1, nspec
       call get_indexed_namelist_unit (unit, "species_parameters", is)
       z = 1
       mass = 1.0
       dens = 1.0
       temp = 1.0
       tprim = -999.9
       fprim = -999.9
       d2ndr2 = 0.0
       d2Tdr2 = 0.0
       type = "default"
       read (unit=unit, nml=species_parameters)
       close (unit=unit)

       spec(is)%z = z
       spec(is)%mass = mass
       spec(is)%dens = dens
       spec(is)%temp = temp
       spec(is)%tprim = tprim
       spec(is)%fprim = fprim
       ! this is (1/n_s)*d^2 n_s / drho^2
       spec(is)%d2ndr2 = d2ndr2
       ! this is (1/T_s)*d^2 T_s / drho^2
       spec(is)%d2Tdr2 = d2Tdr2

       ierr = error_unit()
       call get_option_value (type, typeopts, spec(is)%type, ierr, "type in species_parameters_x")
    end do

  end subroutine read_species_stella

  subroutine broadcast_parameters

    use mp, only: broadcast

    implicit none

    integer :: is

    do is = 1, nspec
       call broadcast (spec(is)%z)
       call broadcast (spec(is)%mass)
       call broadcast (spec(is)%dens)
       call broadcast (spec(is)%temp)
       call broadcast (spec(is)%tprim)
       call broadcast (spec(is)%fprim)
       call broadcast (spec(is)%vnew)
       call broadcast (spec(is)%d2ndr2)
       call broadcast (spec(is)%d2Tdr2)
       call broadcast (spec(is)%type)

       spec(is)%stm = sqrt(spec(is)%temp/spec(is)%mass)
       spec(is)%zstm = spec(is)%z/sqrt(spec(is)%temp*spec(is)%mass)
       spec(is)%tz = spec(is)%temp/spec(is)%z
       spec(is)%zt = spec(is)%z/spec(is)%temp
       spec(is)%smz = abs(sqrt(spec(is)%temp*spec(is)%mass)/spec(is)%z)
    end do

  end subroutine broadcast_parameters

  pure function has_electron_species (spec)
    use common_types, only: spec_type
    implicit none
    type (spec_type), dimension (:), intent (in) :: spec
    logical :: has_electron_species
    has_electron_species = any(spec%type == electron_species)
  end function has_electron_species

  pure function has_slowing_down_species (spec)
    use common_types, only: spec_type
    implicit none
    type (spec_type), dimension (:), intent (in) :: spec
    logical :: has_slowing_down_species
    has_slowing_down_species = any(spec%type == slowing_down_species)
  end function has_slowing_down_species

  subroutine finish_species

    implicit none

    deallocate (spec)

    initialized = .false.

  end subroutine finish_species

!   subroutine reinit_species (ntspec, dens, temp, fprim, tprim, nu)

!     use mp, only: broadcast, proc0

!     implicit none

!     integer, intent (in) :: ntspec
!     real, dimension (:), intent (in) :: dens, fprim, temp, tprim, nu

!     integer :: is
!     logical, save :: first = .true.

!     if (first) then
!        if (nspec == 1) then
!           ions = 1
!           electrons = 0
!           impurity = 0
!        else
!           ! if 2 or more species in GS2 calculation, figure out which is main ion
!           ! and which is electron via mass (main ion mass assumed to be one)
!           do is = 1, nspec
!              if (abs(spec(is)%mass-1.0) <= epsilon(0.0)) then
!                 ions = is
!              else if (spec(is)%mass < 0.3) then
!                 ! for electrons, assuming electrons are at least a factor of 3 less massive
!                 ! than main ion and other ions are no less than 30% the mass of the main ion
!                 electrons = is
!              else if (spec(is)%mass > 1.0 + epsilon(0.0)) then
!                 impurity = is
!              else
!                 if (proc0) write (*,*) &
!                      "Error: TRINITY requires the main ions to have mass 1", &
!                      "and the secondary ions to be impurities (mass > 1)"
!                 stop
!              end if
!           end do
!        end if
!        first = .false.
!     end if

!     if (proc0) then

!        nspec = ntspec

!        ! TRINITY passes in species in following order: main ion, electron, impurity (if present)

!        ! for now, hardwire electron density to be reference density
!        ! main ion temperature is reference temperature
!        ! main ion mass is assumed to be the reference mass
       
!        ! if only 1 species in the GS2 calculation, it is assumed to be main ion
!        ! and ion density = electron density
!        if (nspec == 1) then
!           spec(1)%dens = 1.0
!           spec(1)%temp = 1.0
!           spec(1)%fprim = fprim(1)
!           spec(1)%tprim = tprim(1)
! !          spec(1)%vnewk = nu(1)
!        else
!           spec(ions)%dens = dens(1)/dens(2)
!           spec(ions)%temp = 1.0
!           spec(ions)%fprim = fprim(1)
!           spec(ions)%tprim = tprim(1)
! !          spec(ions)%vnewk = nu(1)

!           spec(electrons)%dens = 1.0
!           spec(electrons)%temp = temp(2)/temp(1)
!           spec(electrons)%fprim = fprim(2)
!           spec(electrons)%tprim = tprim(2)
! !          spec(electrons)%vnewk = nu(2)

!           if (nspec > 2) then
!              spec(impurity)%dens = dens(3)/dens(2)
!              spec(impurity)%temp = temp(3)/temp(1)
!              spec(impurity)%fprim = fprim(3)
!              spec(impurity)%tprim = tprim(3)
! !             spec(impurity)%vnewk = nu(3)
!           end if
!        end if

!        do is = 1, nspec
!           spec(is)%stm = sqrt(spec(is)%temp/spec(is)%mass)
!           spec(is)%zstm = spec(is)%z/sqrt(spec(is)%temp*spec(is)%mass)
!           spec(is)%tz = spec(is)%temp/spec(is)%z
!           spec(is)%zt = spec(is)%z/spec(is)%temp
!           spec(is)%smz = abs(sqrt(spec(is)%temp*spec(is)%mass)/spec(is)%z)

! !          write (*,100) 'reinit_species', rhoc_ms, spec(is)%temp, spec(is)%fprim, &
! !               spec(is)%tprim, spec(is)%vnewk, real(is)
!        end do

!     end if

! !100 format (a15,9(1x,1pg18.11))

!     call broadcast (nspec)

!     do is = 1, nspec
!        call broadcast (spec(is)%dens)
!        call broadcast (spec(is)%temp)
!        call broadcast (spec(is)%fprim)
!        call broadcast (spec(is)%tprim)
! !       call broadcast (spec(is)%vnewk)
!        call broadcast (spec(is)%stm)
!        call broadcast (spec(is)%zstm)
!        call broadcast (spec(is)%tz)
!        call broadcast (spec(is)%zt)
!        call broadcast (spec(is)%smz)
!     end do

!   end subroutine reinit_species

!   subroutine init_trin_species (ntspec_in, dens_in, temp_in, fprim_in, tprim_in, nu_in)

!     implicit none

!     integer, intent (in) :: ntspec_in
!     real, dimension (:), intent (in) :: dens_in, fprim_in, temp_in, tprim_in, nu_in

!     if (.not. allocated(temp_trin)) then
!        allocate (dens_trin(size(dens_in)))
!        allocate (fprim_trin(size(fprim_in)))
!        allocate (temp_trin(size(temp_in)))
!        allocate (tprim_trin(size(tprim_in)))
!        allocate (nu_trin(size(nu_in)))
!     end if

!     ntspec_trin = ntspec_in
!     dens_trin = dens_in
!     temp_trin = temp_in
!     fprim_trin = fprim_in
!     tprim_trin = tprim_in
!     nu_trin = nu_in

!   end subroutine init_trin_species

end module species
