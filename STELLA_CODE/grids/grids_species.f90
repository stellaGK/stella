module grids_species

   use stella_common_types, only: spec_type
   
   ! Read the parameters for <species_option_switch> from namelist_species.f90
   use namelist_species, only: species_option_stella
   use namelist_species, only: species_option_inputprofs
   use namelist_species, only: species_option_euterpe
   use namelist_species, only: species_option_multibox
   
   ! Read the parameters for <type> from namelist_species.f90
   use namelist_species, only: ion_species
   use namelist_species, only: electron_species
   use namelist_species, only: slowing_down_species
   use namelist_species, only: tracer_species

   implicit none

   ! Make routines available to other modules
   public :: init_species, finish_species
   public :: read_parameters_species
   public :: reinit_species
   public :: communicate_species_multibox
   
   ! Although the parameters are available through namelist_species, 
   ! make them available through grids_species as well
   public :: ion_species, electron_species
   public :: slowing_down_species, tracer_species
   public :: species_option_stella, species_option_inputprofs
   public :: species_option_euterpe, species_option_multibox
   
   public :: nspec, spec, pfac
   public :: has_electron_species, has_slowing_down_species
   public :: ions, electrons, impurity
   public :: modified_adiabatic_electrons, adiabatic_electrons
   
   ! Allow other routines to check whether we are using multibox
   public :: species_option_switch

   private
   
   integer :: species_option_switch

   integer :: nspec
   logical :: read_profile_variation, write_profile_variation
   logical :: ecoll_zeff
   logical :: modified_adiabatic_electrons, adiabatic_electrons

   type(spec_type), dimension(:), allocatable :: spec

   integer :: ions, electrons, impurity
   real :: pfac
   
   logical :: initialised_init_species = .false.
   logical :: initialised_read_species = .false.

contains

!###############################################################################
!################################ READ NAMELIST ################################
!###############################################################################

   subroutine read_parameters_species
   
      use mp, only: proc0
      use namelist_species, only: read_namelist_species_options
      use namelist_species, only: read_namelist_species_stella
      use grids_species_from_euterpe, only: read_species_euterpe
      use geometry_inputprofiles_interface, only: read_inputprof_spec

      implicit none
      
      !-------------------------------------------------------------------------
      
      ! Only initialise once
      if (initialised_read_species) return
      initialised_read_species = .true.

      ! Read the "species_options" namelists in the input file
      call read_namelist_species_options(nspec, species_option_switch, &
         read_profile_variation, write_profile_variation, ecoll_zeff)
         
      ! Broadcast the parameters to all processors
      call broadcast_species_options
         
      ! Allocate species now that we know <nspec>
      allocate (spec(nspec))
         
      ! Read the "species_parameters_i" and "euterpe_parameters" namelists in the input file
       call read_namelist_species_stella(nspec, spec)
       
      ! Read the "euterpe_parameters" namelists in the input file or an input profile
      if (proc0) then
         if (species_option_switch == species_option_inputprofs) then
            call read_inputprof_spec(nspec, spec)
         end if
         if (species_option_switch == species_option_euterpe) then
            call read_species_euterpe(nspec, spec)
         end if
      end if
      
   contains
   
      subroutine broadcast_species_options

         use mp, only: broadcast

         implicit none

         call broadcast(nspec)
         call broadcast(read_profile_variation)
         call broadcast(write_profile_variation)
         call broadcast(species_option_switch)
         call broadcast(ecoll_zeff)

      end subroutine broadcast_species_options
      
   end subroutine read_parameters_species

!###############################################################################
!############################### INITIALISE Z GRID #############################
!###############################################################################

   !****************************************************************************
   !                                INITIALISE SPECIES
   !****************************************************************************
   subroutine init_species

      use mp, only: proc0, broadcast, mp_abort
      use parameters_physics, only: vnew_ref, zeff
      use parameters_multibox, only: include_pressure_variation
      use parameters_physics, only: adiabatic_option_switch, adiabatic_option_fieldlineavg

      use parameters_physics, only: full_flux_surface
      implicit none

      integer :: is, is2

      if (initialised_init_species) return
      initialised_init_species = .true.

      !this will be called by the central box in stella.f90 after
      !ktgrids is set up as we need to know the radial box size
      if (proc0) then
         if (species_option_switch == species_option_multibox) then
            call communicate_species_multibox
         end if
      end if

      if (proc0) then
         if (ecoll_zeff) then
            ! AVB: only intra-species collisions, account for e-i and e-impurity collisions using zeff
            do is = 1, nspec
               ! initialize nu_ss' = 0 for all s'
               spec(is)%vnew = 0.
               ! FLAG -- only contains self-collisions at the moment
               spec(is)%vnew(is) = vnew_ref * spec(is)%dens * spec(is)%z**4 &
                                   / (sqrt(spec(is)%mass) * spec(is)%temp**1.5)
               ! include electron-ion collisions
               if (spec(is)%type == electron_species) then
                  spec(is)%vnew(is) = spec(is)%vnew(is) * (1.+zeff)
               end if
            end do
         else
            ! AVB: full intra- and inter-species collision frequencies
            do is = 1, nspec
               do is2 = 1, nspec
                  if (spec(is)%type == electron_species) then
                     spec(is)%vnew(is2) = vnew_ref * spec(is2)%dens * spec(is)%z**2 * spec(is2)%z**2 &
                                          / (sqrt(spec(is)%mass) * spec(is)%temp**1.5)
                  else if ((spec(is)%type == ion_species) .and. (spec(is2)%type == ion_species)) then
                     spec(is)%vnew(is2) = vnew_ref * spec(is2)%dens * spec(is)%z**2 * spec(is2)%z**2 &
                                          / (sqrt(spec(is)%mass) * spec(is)%temp**1.5)
                  else if ((spec(is)%type == ion_species) .and. (spec(is2)%type == electron_species)) then
                     spec(is)%vnew(is2) = vnew_ref * spec(is2)%dens * spec(is)%z**2 * spec(is2)%z**2 &
                                          / (sqrt(spec(is)%mass) * spec(is)%temp**1.5)
                  end if
               end do
            end do
         end if

         call dump_species_input

      end if

      pfac = 1.0
      if (.not. include_pressure_variation) pfac = 0.0

      call broadcast_parameters

      ! set flag adiabatic_electrons to true if no kinetic electron species evolved
      adiabatic_electrons = .not. has_electron_species(spec)
      ! set flag modified_adiabatic_electrons to true if no kinetic electron species evolved
      ! and field-line-avg chosen as the adiabatic option
      modified_adiabatic_electrons = adiabatic_electrons &
                                     .and. adiabatic_option_switch == adiabatic_option_fieldlineavg

      if(has_electron_species(spec) .and. full_flux_surface) then
         write (*,*) 'full_flux_surface is not set up for kinetic electrons yet'
         call mp_abort('full_flux_surface is not set up for kinetic electrons yet')
      end if

   end subroutine init_species

   subroutine broadcast_parameters

      use mp, only: broadcast

      implicit none

      integer :: is

      do is = 1, nspec
         call broadcast(spec(is)%z)
         call broadcast(spec(is)%mass)
         call broadcast(spec(is)%dens)
         call broadcast(spec(is)%temp)
         call broadcast(spec(is)%tprim)
         call broadcast(spec(is)%fprim)
         call broadcast(spec(is)%vnew)
         call broadcast(spec(is)%d2ndr2)
         call broadcast(spec(is)%d2Tdr2)
         call broadcast(spec(is)%dens_psi0)
         call broadcast(spec(is)%temp_psi0)
         call broadcast(spec(is)%bess_fac)
         call broadcast(spec(is)%type)

         spec(is)%stm = sqrt(spec(is)%temp / spec(is)%mass)
         spec(is)%zstm = spec(is)%z / sqrt(spec(is)%temp * spec(is)%mass)
         spec(is)%tz = spec(is)%temp / spec(is)%z
         spec(is)%zt = spec(is)%z / spec(is)%temp
         spec(is)%smz = abs(sqrt(spec(is)%temp * spec(is)%mass) / spec(is)%z)

         spec(is)%stm_psi0 = sqrt(spec(is)%temp_psi0 / spec(is)%mass)
         spec(is)%zstm_psi0 = spec(is)%z / sqrt(spec(is)%temp_psi0 * spec(is)%mass)
         spec(is)%tz_psi0 = spec(is)%temp_psi0 / spec(is)%z
         spec(is)%zt_psi0 = spec(is)%z / spec(is)%temp_psi0
         spec(is)%smz_psi0 = abs(sqrt(spec(is)%temp_psi0 * spec(is)%mass) / spec(is)%z)
      end do

   end subroutine broadcast_parameters

   !****************************************************************************
   !                                FUNCTIONS
   !****************************************************************************
   pure function has_electron_species(spec)
      use stella_common_types, only: spec_type
      implicit none
      type(spec_type), dimension(:), intent(in) :: spec
      logical :: has_electron_species
      has_electron_species = any(spec%type == electron_species)
   end function has_electron_species

   pure function has_slowing_down_species(spec)
      use stella_common_types, only: spec_type
      implicit none
      type(spec_type), dimension(:), intent(in) :: spec
      logical :: has_slowing_down_species
      has_slowing_down_species = any(spec%type == slowing_down_species)
   end function has_slowing_down_species

   !****************************************************************************
   !                                FINISH SPECIES
   !****************************************************************************
   subroutine finish_species

      implicit none

      deallocate (spec)

      initialised_init_species = .false.

   end subroutine finish_species

   subroutine reinit_species(ntspec, dens, temp, fprim, tprim, bess_fac)

      use mp, only: broadcast, proc0, mp_abort

      implicit none

      integer, intent(in) :: ntspec
      real, dimension(:), intent(in) :: dens, fprim, temp, tprim, bess_fac

      integer :: is
      logical, save :: first = .true.

      if (first) then
         if (nspec == 1) then
            ions = 1
            electrons = 0
            impurity = 0
         else
            ! if 2 or more species in GS2 calculation, figure out which is main ion
            ! and which is electron via mass (main ion mass assumed to be one)
            do is = 1, nspec
               if (abs(spec(is)%mass - 1.0) <= epsilon(0.0)) then
                  ions = is
               else if (spec(is)%mass < 0.3) then
                  ! for electrons, assuming electrons are at least a factor of 3 less massive
                  ! than main ion and other ions are no less than 30% the mass of the main ion
                  electrons = is
               else if (spec(is)%mass > 1.0 + epsilon(0.0)) then
                  impurity = is
               else
                  if (proc0) write (*, *) &
                     "Error: TRINITY requires the main ions to have mass 1", &
                     "and the secondary ions to be impurities (mass > 1)"
                  call mp_abort('TRINITY requires the main ions to have mass 1 and the secondary ions mass > 1')
               end if
            end do
         end if
         first = .false.
      end if

      if (proc0) then

         nspec = ntspec

         ! Species are passed in following order: main ion, electron, impurity (if present)
         if (nspec == 1) then
            spec(1)%dens = dens(1)
            spec(1)%temp = temp(1)
            spec(1)%fprim = fprim(1)
            spec(1)%tprim = tprim(1)
            spec(1)%bess_fac = bess_fac(1)
         else
            spec(ions)%dens = dens(1)
            spec(ions)%temp = temp(1)
            spec(ions)%fprim = fprim(1)
            spec(ions)%tprim = tprim(1)
            spec(ions)%bess_fac = bess_fac(1)

            spec(electrons)%dens = dens(2)
            spec(electrons)%temp = temp(2)
            spec(electrons)%fprim = fprim(2)
            spec(electrons)%tprim = tprim(2)
            spec(electrons)%bess_fac = bess_fac(2)

            if (nspec > 2) then
               spec(impurity)%dens = dens(3)
               spec(impurity)%temp = temp(3)
               spec(impurity)%fprim = fprim(3)
               spec(impurity)%tprim = tprim(3)
               spec(impurity)%bess_fac = bess_fac(3)
            end if
         end if

         do is = 1, nspec
            spec(is)%stm = sqrt(spec(is)%temp / spec(is)%mass)
            spec(is)%zstm = spec(is)%z / sqrt(spec(is)%temp * spec(is)%mass)
            spec(is)%tz = spec(is)%temp / spec(is)%z
            spec(is)%zt = spec(is)%z / spec(is)%temp
            spec(is)%smz = abs(sqrt(spec(is)%temp * spec(is)%mass) / spec(is)%z)

            !          write (*,100) 'reinit_species', rhoc_ms, spec(is)%temp, spec(is)%fprim, &
            !               spec(is)%tprim, spec(is)%vnewk, real(is)
         end do

         call dump_species_input

      end if

      !100 format (a15,9(1x,1pg18.11))

      call broadcast(nspec)

      do is = 1, nspec
         call broadcast(spec(is)%dens)
         call broadcast(spec(is)%temp)
         call broadcast(spec(is)%fprim)
         call broadcast(spec(is)%tprim)
         call broadcast(spec(is)%bess_fac)
         call broadcast(spec(is)%stm)
         call broadcast(spec(is)%zstm)
         call broadcast(spec(is)%tz)
         call broadcast(spec(is)%zt)
         call broadcast(spec(is)%smz)
      end do

   end subroutine reinit_species

   subroutine communicate_species_multibox(dr_m, dr_p)
      use job_manage, only: njobs
      use mp, only: job, scope, mp_abort, &
                    crossdomprocs, subprocs, &
                    send, receive

      implicit none

      real, optional, intent(in) :: dr_m, dr_p

      real, dimension(:), allocatable :: dens, ldens, ltemp, lfprim, ltprim
      real, dimension(:), allocatable :: temp, rdens, rtemp, rfprim, rtprim

      integer :: i

      allocate (dens(nspec))
      allocate (temp(nspec))
      allocate (ldens(nspec))
      allocate (ltemp(nspec))
      allocate (lfprim(nspec))
      allocate (ltprim(nspec))
      allocate (rdens(nspec))
      allocate (rtemp(nspec))
      allocate (rfprim(nspec))
      allocate (rtprim(nspec))

      if (job == 1) then

         ! recall that fprim and tprim are the negative gradients
         ldens = spec%dens * (1.0 - dr_m * spec%fprim)! + 0.5*dr_m**2*spec%d2ndr2)
         ltemp = spec%temp * (1.0 - dr_m * spec%tprim)! + 0.5*dr_m**2*spec%d2Tdr2)
         lfprim = (spec%fprim - dr_m * spec%d2ndr2) * (spec%dens / ldens)
         ltprim = (spec%tprim - dr_m * spec%d2Tdr2) * (spec%temp / ltemp)

         rdens = spec%dens * (1.0 - dr_p * spec%fprim)! + 0.5*dr_p**2*spec%d2ndr2)
         rtemp = spec%temp * (1.0 - dr_p * spec%tprim)! + 0.5*dr_p**2*spec%d2Tdr2)
         rfprim = (spec%fprim - dr_p * spec%d2ndr2) * (spec%dens / rdens)
         rtprim = (spec%tprim - dr_p * spec%d2Tdr2) * (spec%temp / rtemp)

         do i = 1, nspec
            if (ldens(i) < 0 .or. ltemp(i) < 0 .or. &
                rdens(i) < 0 .or. rtemp(i) < 0) then
               call mp_abort('Negative n/T encountered. Try reducing rhostar.')
            end if
         end do
      end if

      call scope(crossdomprocs)

      if (job == 1) then
         call send(ldens, 0, 120)
         call send(ltemp, 0, 121)
         call send(lfprim, 0, 122)
         call send(ltprim, 0, 123)
         call send(spec%dens, 0, 124)
         call send(spec%temp, 0, 125)
         call send(rdens, njobs - 1, 130)
         call send(rtemp, njobs - 1, 131)
         call send(rfprim, njobs - 1, 132)
         call send(rtprim, njobs - 1, 133)
         call send(spec%dens, njobs - 1, 134)
         call send(spec%temp, njobs - 1, 135)
      elseif (job == 0) then
         call receive(ldens, 1, 120)
         call receive(ltemp, 1, 121)
         call receive(lfprim, 1, 122)
         call receive(ltprim, 1, 123)
         call receive(dens, 1, 124)
         call receive(temp, 1, 125)
         spec%dens = ldens
         spec%temp = ltemp
         spec%fprim = lfprim
         spec%tprim = ltprim
         spec%dens_psi0 = dens
         spec%temp_psi0 = temp
      elseif (job == njobs - 1) then
         call receive(rdens, 1, 130)
         call receive(rtemp, 1, 131)
         call receive(rfprim, 1, 132)
         call receive(rtprim, 1, 133)
         call receive(dens, 1, 134)
         call receive(temp, 1, 135)
         spec%dens = rdens
         spec%temp = rtemp
         spec%fprim = rfprim
         spec%tprim = rtprim
         spec%dens_psi0 = dens
         spec%temp_psi0 = temp
      end if

      call scope(subprocs)

      deallocate (dens)
      deallocate (temp)
      deallocate (ldens)
      deallocate (ltemp)
      deallocate (lfprim)
      deallocate (ltprim)
      deallocate (rdens)
      deallocate (rtemp)
      deallocate (rfprim)
      deallocate (rtprim)

   end subroutine communicate_species_multibox

   subroutine dump_species_input

      use file_utils, only: run_name

      implicit none

      integer :: is
      character(300) :: filename

      filename = trim(trim(run_name)//'.species.input')
      open (1003, file=filename, status='unknown')
      write (1003, '(9a12,a9)') '#1.z', '2.mass', '3.dens', &
         '4.temp', '5.tprim', '6.fprim', '7.vnewss', &
         '8.dens_psi0', '9.temp_psi0', '11.type'
      do is = 1, nspec
         write (1003, '(9e12.4,i9)') spec(is)%z, spec(is)%mass, &
            spec(is)%dens, spec(is)%temp, spec(is)%tprim, &
            spec(is)%fprim, spec(is)%vnew(is), &
            spec(is)%dens_psi0, spec(is)%temp_psi0, &
            spec(is)%type
      end do
      close (1003)

   end subroutine dump_species_input

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

end module grids_species
