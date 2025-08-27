!###############################################################################
!                                MODULE SPECIES
!###############################################################################
! This module reads the number of species and the specifics of each specie.
!###############################################################################
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
   public :: communicate_species_multibox
   
   ! Although the parameters are available through namelist_species, 
   ! make them available through grids_species as well
   public :: ion_species, electron_species
   public :: slowing_down_species, tracer_species
   public :: species_option_stella, species_option_inputprofs
   public :: species_option_euterpe, species_option_multibox
   
   ! Make the following parameters public
   public :: nspec, spec
   public :: has_electron_species, has_slowing_down_species
   public :: modified_adiabatic_electrons, adiabatic_electrons
   
   ! Pressure variation
   public :: pfac
   
   ! Allow other routines to check whether we are using multibox
   public :: species_option_switch

   private
   
   integer :: species_option_switch

   integer :: nspec
   logical :: read_profile_variation, write_profile_variation
   logical :: modified_adiabatic_electrons, adiabatic_electrons

   type(spec_type), dimension(:), allocatable :: spec
   
   ! Pressure variation
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
         read_profile_variation, write_profile_variation)
         
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

      end subroutine broadcast_species_options
      
   end subroutine read_parameters_species

!###############################################################################
!########################### INITIALISE SPECIES GRID ###########################
!###############################################################################

   !****************************************************************************
   !                                INITIALISE SPECIES
   !****************************************************************************
   subroutine init_species(ecoll_zeff)

      use mp, only: proc0, broadcast, mp_abort
      use parameters_multibox, only: include_pressure_variation
      use parameters_physics, only: adiabatic_option_switch, adiabatic_option_fieldlineavg
      use parameters_physics, only: full_flux_surface
      
      implicit none

      ! Mote that the dissipation module depends on the species module so we 
      ! can not import <ecoll_zeff> and have to pass it through instead
      logical, intent (in) :: ecoll_zeff

      ! Local variables
      integer :: is, is2
      
      !-------------------------------------------------------------------------

      ! Only initialise once
      if (initialised_init_species) return
      initialised_init_species = .true.
      
      ! Based on <ecoll_zeff> set the collision frequencies
      call init_collision_frequencies

      ! Pressure variation
      pfac = 1.0
      if (.not. include_pressure_variation) pfac = 0.0

      ! Broadcast the species parameters
      call broadcast_species_parameters

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
      
   contains
     
      !-------------------- Initialise collision frequencies -------------------
      subroutine init_collision_frequencies
      
         ! Parallelisation
         use mp, only: proc0
          
         ! Get the reference collisions frequency
         ! and the effective charge of the plasma
         use parameters_physics, only: vnew_ref
         use parameters_physics, only: zeff

         implicit none
         
         !-------------------------------------------------------------------------
 
         ! Calculate collision frequencies only on the first processor
         if (.not. proc0) return

         ! If <ecoll_zeff> = True, we only consider intra-species collisions, so
         ! account for e-i and e-impurity collisions using zeff.
         if (ecoll_zeff) then
         
            ! Iterate over the species
            do is = 1, nspec
            
               ! Initialise nu_ss' = 0 for all s'
               spec(is)%vnew = 0.
               
               ! Set collision frequency for self collisions, i.e., is=is
               ! FLAG -- only contains self-collisions at the moment
               spec(is)%vnew(is) = vnew_ref * spec(is)%dens * spec(is)%z**4 &
                   / (sqrt(spec(is)%mass) * spec(is)%temp**1.5)
                   
               ! Set collision frequency for electron-ion collisions
               ! FLAG -- aren't these electron-electron collisions?
               if (spec(is)%type == electron_species) then
                  spec(is)%vnew(is) = spec(is)%vnew(is) * (1.+zeff)
               end if
               
            end do
            
         ! If <ecoll_zeff> = False, we consider full intra- and inter-species collision frequencies
         else
         
            ! Iterate over the species
            do is = 1, nspec
               do is2 = 1, nspec
               
                  ! Set electron-electron collision frequency
                  if (spec(is)%type == electron_species) then
                     spec(is)%vnew(is2) = vnew_ref * spec(is2)%dens * spec(is)%z**2 * spec(is2)%z**2 &
                         / (sqrt(spec(is)%mass) * spec(is)%temp**1.5)
                         
                  ! Set ion-ion collision frequency
                  else if ((spec(is)%type == ion_species) .and. (spec(is2)%type == ion_species)) then
                     spec(is)%vnew(is2) = vnew_ref * spec(is2)%dens * spec(is)%z**2 * spec(is2)%z**2 &
                         / (sqrt(spec(is)%mass) * spec(is)%temp**1.5)
                         
                  ! Set ion-electron collision frequency
                  else if ((spec(is)%type == ion_species) .and. (spec(is2)%type == electron_species)) then
                     spec(is)%vnew(is2) = vnew_ref * spec(is2)%dens * spec(is)%z**2 * spec(is2)%z**2 &
                         / (sqrt(spec(is)%mass) * spec(is)%temp**1.5)
                  end if
                  
               end do
            end do
            
         end if

         call write_species_parameters_to_txt_file
      
      end subroutine init_collision_frequencies

      !--------------------------- Broadcast parameters----------------------------
      subroutine broadcast_species_parameters

         use mp, only: broadcast

         implicit none

         integer :: is
         
         !-------------------------------------------------------------------------

         ! Iterate over the species
         do is = 1, nspec
         
            ! Broadcast data from first processor to all processors
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

            ! Calculate extra useful quantities on all  processors
            spec(is)%stm = sqrt(spec(is)%temp / spec(is)%mass)
            spec(is)%zstm = spec(is)%z / sqrt(spec(is)%temp * spec(is)%mass)
            spec(is)%tz = spec(is)%temp / spec(is)%z
            spec(is)%zt = spec(is)%z / spec(is)%temp
            spec(is)%smz = abs(sqrt(spec(is)%temp * spec(is)%mass) / spec(is)%z)

            ! Perform the same calculations on <psi0> for radial variation runs
            spec(is)%stm_psi0 = sqrt(spec(is)%temp_psi0 / spec(is)%mass)
            spec(is)%zstm_psi0 = spec(is)%z / sqrt(spec(is)%temp_psi0 * spec(is)%mass)
            spec(is)%tz_psi0 = spec(is)%temp_psi0 / spec(is)%z
            spec(is)%zt_psi0 = spec(is)%z / spec(is)%temp_psi0
            spec(is)%smz_psi0 = abs(sqrt(spec(is)%temp_psi0 * spec(is)%mass) / spec(is)%z)
            
         end do

      end subroutine broadcast_species_parameters

   end subroutine init_species


!###############################################################################
!################################## FUNCTIONS ##################################
!###############################################################################

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
   
!###############################################################################
!########################## WRITE SPECIES TO TXT FILE ##########################
!###############################################################################
   subroutine write_species_parameters_to_txt_file

      use file_utils, only: run_name

      implicit none

      integer :: is
      character(300) :: filename
         
      !-------------------------------------------------------------------------

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

   end subroutine write_species_parameters_to_txt_file

!###############################################################################
!################################ FINISH SPECIES ###############################
!###############################################################################
   subroutine finish_species

      implicit none

      deallocate (spec)

      initialised_init_species = .false.

   end subroutine finish_species

!###############################################################################
!######################### MULTIBOX - RADIAL VARIATION #########################
!###############################################################################
   subroutine communicate_species_multibox(dr_m, dr_p)
   
      ! Parallelisation
      use job_manage, only: njobs
      use mp, only: job, scope, mp_abort
      use mp, only: crossdomprocs, subprocs
      use mp, only: send, receive

      implicit none

      real, optional, intent(in) :: dr_m, dr_p

      real, dimension(:), allocatable :: dens, ldens, ltemp, lfprim, ltprim
      real, dimension(:), allocatable :: temp, rdens, rtemp, rfprim, rtprim

      integer :: i
         
      !-------------------------------------------------------------------------

      ! Allocate species arrays
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

      ! Iterate over the jobs
      if (job == 1) then

         ! Recall that fprim and tprim are the negative gradients
         ldens = spec%dens * (1.0 - dr_m * spec%fprim)! + 0.5*dr_m**2*spec%d2ndr2)
         ltemp = spec%temp * (1.0 - dr_m * spec%tprim)! + 0.5*dr_m**2*spec%d2Tdr2)
         lfprim = (spec%fprim - dr_m * spec%d2ndr2) * (spec%dens / ldens)
         ltprim = (spec%tprim - dr_m * spec%d2Tdr2) * (spec%temp / ltemp) 
         rdens = spec%dens * (1.0 - dr_p * spec%fprim)! + 0.5*dr_p**2*spec%d2ndr2)
         rtemp = spec%temp * (1.0 - dr_p * spec%tprim)! + 0.5*dr_p**2*spec%d2Tdr2)
         rfprim = (spec%fprim - dr_p * spec%d2ndr2) * (spec%dens / rdens)
         rtprim = (spec%tprim - dr_p * spec%d2Tdr2) * (spec%temp / rtemp)

         ! Iterate over the species
         do i = 1, nspec
            if (ldens(i) < 0 .or. ltemp(i) < 0 .or. &
                rdens(i) < 0 .or. rtemp(i) < 0) then
               call mp_abort('Negative n/T encountered. Try reducing rhostar.')
            end if
         end do
         
      end if

      ! Parallelisation
      call scope(crossdomprocs)

      ! Parallelisation
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

      ! Parallelisation
      call scope(subprocs)

      ! Deallocate arrays
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

end module grids_species
