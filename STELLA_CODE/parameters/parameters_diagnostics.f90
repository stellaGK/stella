module parameters_diagnostics

   use debug_flags, only: debug => diagnostics_parameters

   implicit none

   ! Routine to read "diagnostics_knobs" in the input file
   public :: read_diagnostics_knobs 

   ! Variables used to write diagnostics
   public :: nwrite, nsave, navg, nc_mult
   public :: save_for_restart, autostop, write_all
   
   ! Write time traces
   public :: write_phi2_vs_time
   public :: write_apar2_vs_time
   public :: write_bpar2_vs_time
   public :: write_fluxes_vs_time

   ! Write potential in <diagnostics_potential>
   public :: write_phi_vs_kxkyz
   public :: write_phi2_vs_kxky

   ! Write omega in <diagnostics_omega>
   public :: write_omega_vs_kxky
   public :: write_omega_avg_vs_kxky

   ! Write fluxes in <diagnostics_fluxes>
   public :: write_fluxes_kxkyz
   public :: write_fluxes_kxky
   public :: write_radial_fluxes
   public :: flux_norm 

   ! Write distribution in <diagnostics_distribution>
   public :: write_g2_vs_vpamus
   public :: write_g2_vs_zvpas
   public :: write_g2_vs_zmus 
   public :: write_g2_vs_kxkyzs 
   public :: write_g2_vs_zvpamus 
   public :: write_distribution_g
   public :: write_distribution_h
   public :: write_distribution_f

   ! Write moments in <diagnostics_moments>
   public :: write_radial_moments
   public :: write_moments  

   private

   ! Variables used to write diagnostics
   integer :: nwrite, nsave, navg, nc_mult
   logical :: save_for_restart, autostop, write_all
   
   ! Write time traces
   logical :: write_phi2_vs_time
   logical :: write_apar2_vs_time
   logical :: write_bpar2_vs_time
   logical :: write_fluxes_vs_time

   ! Write potential in <diagnostics_potential>
   logical :: write_phi_vs_kxkyz
   logical :: write_phi2_vs_kxky

   ! Write omega in <diagnostics_omega>
   logical :: write_omega_vs_kxky
   logical :: write_omega_avg_vs_kxky

   ! Write fluxes in <diagnostics_fluxes>
   logical :: write_fluxes_kxkyz
   logical :: write_fluxes_kxky
   logical :: write_radial_fluxes
   logical :: flux_norm 

   ! Write distribution in <diagnostics_distribution>
   logical :: write_g2_vs_vpamus
   logical :: write_g2_vs_zvpas
   logical :: write_g2_vs_zmus 
   logical :: write_g2_vs_kxkyzs 
   logical :: write_g2_vs_zvpamus 
   logical :: write_distribution_g
   logical :: write_distribution_h
   logical :: write_distribution_f

   ! Write moments in <diagnostics_moments>
   logical :: write_radial_moments
   logical :: write_moments     

contains

   !============================================================================
   !====================== READ AND BROADCAST INPUT FILE =======================
   !============================================================================ 
   ! Read-in the parameters from the namelist "diagnostics_knobs" in the 
   ! input file and broadcast the parameters to all the processors.
   ! Gets called in the <init_stella> subroutine in the <stella> module. 
   subroutine read_diagnostics_knobs
 
      use mp, only: proc0

      implicit none
         
      ! Logical old variables for backwards compatibility
      logical :: write_phi_vs_time, write_apar_vs_time, write_bpar_vs_time 
      logical :: write_kspectra, write_gvmus, write_gzvs, write_omega
      
      ! Set the default parameters, read the namelist "diagnostics_knobs" 
      ! in the input file and broadcast the parameters to all processors
      if (proc0) call set_default_parameters 
      if (proc0) call read_input_file 
      call broadcast_parameters 
      
   contains 
   
   
      !**********************************************************************
      !                        SET DEFAULT PARAMETERS                       !
      !********************************************************************** 
      ! Write the time traces by default, while we will not write any higher 
      ! dimensional data by default to have a small netCDF file.
      !********************************************************************** 
      subroutine set_default_parameters

         use parameters_physics, only: radial_variation
         use parameters_physics, only: nonlinear

         implicit none
         
         ! Track code 
         if (debug) write (*, *) 'read_diagnostics_parameters::set_default_parameters'

         ! Write data to the ascii files at every <nwrite> time steps.
         ! Write data to the netcdf file every <nwrite*nc_mult> time steps.
         nwrite = 50 
         nc_mult = 1

         ! Save gvmu(vpa, nmu, i[s,kx,ky,z]) at every <nsave> time steps so that we can restart the simulation.
         ! From <gvmu> we can calculate <phi> through quasi-neutrality, so <gvmu> is all that we need to save.
         save_for_restart = .false.
         nsave = -1
         
         ! Togggle all write statements together 
         write_all = .false.
         
         !------------------------------
         !          Potential          !
         !------------------------------ 
         
         ! Write phi2(t), apar2(t), bpar2(t), phi(kx,ky,z), phi2(kx,ky)
         write_phi2_vs_time = .true.
         write_apar2_vs_time = .true.
         write_bpar2_vs_time = .true. 
         write_phi_vs_kxkyz = .false.
         write_phi2_vs_kxky = .false.
          
         !------------------------------
         !    Distribution function    !
         !------------------------------ 
         
         ! Write the distribution functions g, h and f 
         ! Note that these arrays are very large 
         write_g2_vs_vpamus = .false.
         write_g2_vs_zvpas = .false.
         write_g2_vs_zmus = .false.
         write_g2_vs_kxkyzs = .false.
         write_g2_vs_zvpamus = .false.
         write_distribution_g = .true.
         write_distribution_h = .false.
         write_distribution_f = .false. 
         
         !------------------------------
         !            Omega            !
         !------------------------------ 
         
         ! For linear simulations write omega by default 
         write_omega_vs_kxky = .not. nonlinear
         
         ! We calculate running averages for omega over <navg> time points
         write_omega_avg_vs_kxky = .false.
         navg = 50

         ! Stop linear simulations when gamma is constant (careful since we won't catch jumpers!)
         ! It will check the growth rate gamma over <navg> time steps, hence we need 
         ! <omega_vs_kykx> over the last <navg> time steps, written by <write_omega_avg_vs_kxky>
         autostop = .false.
                  
         !------------------------------
         !           Fluxes            !
         !------------------------------ 
         
         ! Write the particle flux, heat flux and momentum flux
         write_fluxes_vs_time = .true. 
         write_fluxes_kxkyz = .false.   
         write_fluxes_kxky = .false.   

         ! Flux definition with an extra factor 1/<nabla rho> in front.
         ! Toggle to include or not include the factor 1/<∇̃ρ>_ψ in the flux definition
         flux_norm = .true.
         
         !------------------------------
         !           Moments           !
         !------------------------------ 
         
         ! Write the density, temperature and upar
         write_moments = .false.
         
         !------------------------------
         !      Radial variation       !
         !------------------------------ 
         
         ! If <radial_variation> = True, automatically write the corresponding data
         write_radial_fluxes = radial_variation
         write_radial_moments = radial_variation
         
         !------------------------------
         !   Backwards compatibility   !
         !------------------------------ 
         
         ! Backwards compatibility, these flags have been removed
         write_phi_vs_time = .false.
         write_apar_vs_time = .false.
         write_bpar_vs_time = .false.
         write_gvmus = .false.
         write_gzvs = .false. 
         write_kspectra = .false.   
         write_omega = .false.

      end subroutine set_default_parameters 
   
      !**********************************************************************
      !                           READ INPUT FILE                           !
      !**********************************************************************
      ! Define default parameters for the <diagnostics> module and overwrite them
      ! with those defined in the namelist "diagnostics_knobs" in the input file.
      ! Gets called in the <read_diagnostics_knobs> subroutine above. 
      !**********************************************************************
            
      subroutine read_input_file

         use file_utils, only: input_unit_exist
         use parameters_physics, only: nonlinear

         implicit none

         logical :: exist
         integer :: in_file 

         ! Define the namelist "diagnostics_knobs" in the input file.
         ! TODO-HT Change <stella_diagnostics_knobs> to <diagnostics> (Change all input files in tests)
         ! and mp_abort if <stella_diagnostics_knobs> is in the input file
         namelist /stella_diagnostics_knobs/ nwrite, navg, nsave, autostop, save_for_restart, flux_norm, nc_mult, &
            write_phi2_vs_time, write_apar2_vs_time, write_bpar2_vs_time, write_fluxes_vs_time, &
            write_phi_vs_kxkyz, write_g2_vs_vpamus, write_g2_vs_zvpas, write_g2_vs_zmus, &
            write_g2_vs_kxkyzs, write_g2_vs_zvpamus, write_distribution_g, write_distribution_h, write_distribution_f, &
            write_omega_vs_kxky, write_omega_avg_vs_kxky, write_phi2_vs_kxky, write_moments, write_radial_fluxes, &
            write_radial_moments, write_fluxes_kxkyz, write_fluxes_kxky, write_all, flux_norm, nc_mult, &
            ! Backwards compatibility for old stella code
            write_omega, write_phi_vs_time, write_apar_vs_time, write_bpar_vs_time, &
            write_kspectra, write_gvmus, write_gzvs
            
         !-------------------------------------------------------------------
         
         ! Track code 
         if (debug) write (*, *) 'read_diagnostics_parameters::read_input_file'
            
         ! Read the namelist "stella_diagnostics_knobs" in the input file and overwrite the default variables
         in_file = input_unit_exist("stella_diagnostics_knobs", exist)
         if (exist) read (unit=in_file, nml=stella_diagnostics_knobs)

         ! If <save_for_restart> = False then we need <nsave> = -1
         if (.not. save_for_restart) nsave = -1

         ! For nonlinear simulations, don't stop automatically 
         ! and if we want to <autostop> we need the <omega_vs_tkykx> 
         ! which contains <omega_vs_kykx> over the last <navg> time steps
         if (nonlinear) autostop = .false.
         if (autostop) write_omega_avg_vs_kxky = .true.
         
         ! Write all diagnostics
         if (write_all) then 
            write_phi2_vs_time = .true.
            write_apar2_vs_time = .true.
            write_bpar2_vs_time = .true.
            write_fluxes_vs_time = .true. 
            write_omega_vs_kxky = .true.
            write_omega_avg_vs_kxky = .true.
            write_phi_vs_kxkyz = .true.
            write_phi2_vs_kxky = .true.
            write_moments = .true.
            write_fluxes_kxkyz = .true.   
            write_fluxes_kxky = .true.
            write_g2_vs_vpamus = .true.
            write_g2_vs_zvpas = .true.
            write_g2_vs_zmus = .true.
            write_g2_vs_kxkyzs = .true.
            write_g2_vs_zvpamus = .true.
            write_distribution_g = .true.
            write_distribution_h = .true.
            write_distribution_f = .true.  
         end if
         
         ! Backwards compatibility, these flags have been removed
         ! TODO-HT write warnings for backwards compatibility 
         if (write_phi_vs_time) write_phi2_vs_time = .true.
         if (write_apar_vs_time) write_apar2_vs_time = .true.
         if (write_bpar_vs_time) write_bpar2_vs_time = .true.
         if (write_gvmus) write_distribution_g = .true.
         if (write_gvmus) write_g2_vs_vpamus = .true.
         if (write_gzvs) write_distribution_g = .true.
         if (write_gzvs) write_g2_vs_zvpas = .true.
         if (write_kspectra) write_fluxes_kxky = .true. 
         if (write_omega) write_omega_vs_kxky = .true.
         if (write_omega) write_omega_avg_vs_kxky = .true.

      end subroutine read_input_file 
   
      !**********************************************************************
      !                         BROADCAST PARAMETERS                        !
      !**********************************************************************
      subroutine broadcast_parameters
      
         use mp, only: broadcast
      
         implicit none
         
         ! Track code 
         if (debug) write (*, *) 'read_diagnostics_parameters::broadcast_parameters'
       
         ! Writing options 
         call broadcast(nwrite)
         call broadcast(navg)
         call broadcast(nsave)
         call broadcast(nc_mult)
         call broadcast(autostop) 
         call broadcast(save_for_restart)
         call broadcast(flux_norm)
         
         ! Time traces 
         call broadcast(write_phi2_vs_time)
         call broadcast(write_apar2_vs_time)
         call broadcast(write_bpar2_vs_time)
         call broadcast(write_fluxes_vs_time)

         ! Broadcast the variables to all processors
         call broadcast(write_all)
         call broadcast(write_omega_vs_kxky)
         call broadcast(write_omega_avg_vs_kxky)
         call broadcast(write_phi2_vs_kxky)
         call broadcast(write_moments)
         call broadcast(write_phi_vs_kxkyz)
         call broadcast(write_g2_vs_vpamus)
         call broadcast(write_g2_vs_zvpas)
         call broadcast(write_g2_vs_zmus)
         call broadcast(write_g2_vs_kxkyzs)
         call broadcast(write_g2_vs_zvpamus)
         call broadcast(write_distribution_g)
         call broadcast(write_distribution_f)
         call broadcast(write_distribution_h)
         call broadcast(write_radial_fluxes)
         call broadcast(write_radial_moments)
         call broadcast(write_fluxes_kxkyz)
         call broadcast(write_fluxes_kxky)
         
      end subroutine broadcast_parameters

   end subroutine read_diagnostics_knobs

end module parameters_diagnostics
