!###############################################################################
!########################### RANDOM NUMBER GENERATOR ###########################
!###############################################################################
! 
! The random number generator (rng) is used from an external script located at:
!     STELLA_CODE/utils/ran.fpp
! 
! Since the scripts in the utils folder are independent of stella, we use this
! interface to read the <rng_seed> from the stella input file and to initialise
! the random number generator in ran.fpp.
! 
! Later, other modules can directly access the random number generator as,
!     use ran, only: ranf
! 
! The random number generator is used to:
!     - initialise the distribution function with noise if <initialise_distribution_option> = 'noise'
!     - initialise <phase_shift_angle> with noise if <randomize_phase_shift> = True
! 
! Note that if <rng_seed> < 0, the current time will be used as the seed,
! creating random noise. If <rng_seed> > 0, the noise will always be the same.
! For numerical tests, we therefore generally set rng_seed = 1, to ensure that
! the noise initialisation is the same in each simulation.
! 
! The <rng_seed> is read from the "initialise_distribution_noise" namelist.
! 
!###############################################################################

module interface_random_number_generator

   implicit none
   
   ! Make this routine available to other modules
   public :: init_random_number_generator
   
   private
   
contains

   !----------------------------------------------------------------------------
   !-------------------- Initialise random number generator --------------------
   !----------------------------------------------------------------------------
   ! Set up the random number generator (rng) used to initialise the distribution
   ! function g(kx,ky,z,mu,vpa,species) with noise if <initialise_distribution_option> = 'noise', 
   ! and for the <phase_shift_angle> in grids_kxky if <randomize_phase_shift> = True
   ! The <rng_seed> is read from the "initialise_distribution_noise" namelist.
   ! 
   ! Note that if <rng_seed> < 0, the current time will be used as the seed,
   ! creating random noise. If <rng_seed> > 0, the noise will always be the same.
   !----------------------------------------------------------------------------
   subroutine init_random_number_generator
      
      use mp, only: broadcast
         
      implicit none
      
      ! Initialise the rng_seed
      integer :: rng_seed = -1
      
      !----------------------------------------------------------------------
      
      ! Read the <rng_seed> from the input file
      call read_rng_seed_from_initialise_distribution_noise(rng_seed)
      
      ! Broadcast it's value to all processors
      call broadcast(rng_seed)
      
      ! Set up the random number generator on all processors
      call set_up_ranf_module(rng_seed)
      
   contains
   
      !-------------------- Read <rng_seed> from input file --------------------
      subroutine read_rng_seed_from_initialise_distribution_noise(rng_seed)
      
         use mp, only: proc0
         use file_utils, only: input_unit_exist

         implicit none
         
         ! We want to read the rng_seed from the input file 
         integer, intent (in out)  :: rng_seed
         
         ! We need to declare the other variables inside the "initialise_distribution_noise" n
         ! namelist, otherwise, the reading of the namelist will fail
         real :: zf_init = 1.0
         logical :: left = .false.
         logical :: chop_side = .true.
         
         ! Variables needed to read the input file 
         integer :: in_file
         logical :: dexist
         
         ! Variables in the <initialise_distribution_noise> namelist
         namelist /initialise_distribution_noise/ zf_init, left, chop_side, rng_seed
            
         !-------------------------------------------------------------------------

         ! Only read input file on the first processor
         if (.not. proc0) return
         
         ! Overwrite the default input parameters by those specified in the input file
         ! Note that we do not care about <zf_init>, <left> or <chop_side>.
         ! We only want to read in <rng_seed>
         in_file = input_unit_exist('initialise_distribution_noise', dexist)
         if (dexist) read (unit=in_file, nml=initialise_distribution_noise) 
      
      end subroutine read_rng_seed_from_initialise_distribution_noise
      
      !--------------------------- Set up ranf module --------------------------
      subroutine set_up_ranf_module(rng_seed)
      
         use ran, only: get_rnd_seed_length
         use ran, only: init_ranf
         use mp, only: job
         
         implicit none
         
         ! The rng_seed defined in the input file
         integer, intent (in)  :: rng_seed
         
         ! Local variables
         integer, dimension(:), allocatable  :: seed
         integer :: i, n
         
         !----------------------------------------------------------------------
      
         n = get_rnd_seed_length()
         allocate (seed(n))
         if (rng_seed < 0) then
            call init_ranf(.true., seed, job + 2)
         else
            seed = rng_seed + 37 * (/(i - 1, i=1, n)/)
            call init_ranf(.false., seed, job + 2)
         end if
         deallocate (seed)
      
      end subroutine set_up_ranf_module
         
   end subroutine init_random_number_generator
   
end module interface_random_number_generator
