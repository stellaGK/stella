module parameters_diagnostics

   use debug_flags, only: debug => diagnostics_parameters

   implicit none

   ! Routine to read "diagnostics_knobs" in the input file
   public :: read_parameters_diagnostics 

   ! Variables used to write diagnostics
   public :: nwrite, nsave, navg, nc_mult
   public :: save_for_restart, write_all
   
   ! Write time traces
   public :: write_phi2_vs_time
   public :: write_apar2_vs_time
   public :: write_bpar2_vs_time
   public :: write_fluxes_vs_time

   ! Write potential in <diagnostics_potential>
   public :: write_phi_vs_kxkyz
   public :: write_apar_vs_kxkyz
   public :: write_bpar_vs_kxkyz
   public :: write_phi2_vs_kxky
   public :: write_apar2_vs_kxky
   public :: write_bpar2_vs_kxky

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
   logical :: save_for_restart, write_all

   ! Flags to toggle groups 
   logical :: write_all_time_traces
   logical :: write_all_spectra_kxkyz
   logical :: write_all_spectra_kxky
   logical :: write_all_velocity_space
   logical :: write_all_potential
   logical :: write_all_omega
   logical :: write_all_distribution
   logical :: write_all_fluxes
   logical :: write_all_moments
   logical :: write_all_potential_time_traces
   logical :: write_all_potential_spectra
   
   ! Write time traces
   logical :: write_phi2_vs_time
   logical :: write_apar2_vs_time
   logical :: write_bpar2_vs_time
   logical :: write_fluxes_vs_time

   ! Write potential in <diagnostics_potential>
   logical :: write_phi_vs_kxkyz
   logical :: write_apar_vs_kxkyz
   logical :: write_bpar_vs_kxkyz
   logical :: write_phi2_vs_kxky
   logical :: write_apar2_vs_kxky
   logical :: write_bpar2_vs_kxky

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
   subroutine read_parameters_diagnostics
 
      use mp, only: proc0

      implicit none
      
      ! Set the default parameters, read the namelist "diagnostics_knobs" 
      ! in the input file and broadcast the parameters to all processors
      if (proc0) call read_input_file 
      call broadcast_parameters 
      
   contains 
   
      !**********************************************************************
      !                           READ INPUT FILE                           !
      !**********************************************************************
      ! Define default parameters for the <diagnostics> module and overwrite them
      ! with those defined in the namelist "diagnostics_knobs" in the input file.
      ! Gets called in the <read_parameters_diagnostics> subroutine above. 
      !**********************************************************************
            
      subroutine read_input_file

         use namelist_diagnostics, only: &
            read_namelist_diagnostics, read_namelist_diagnostics_potential, &
            read_namelist_diagnostics_omega, read_namelist_diagnostics_distribution, &
            read_namelist_diagnostics_fluxes, read_namelist_diagnostics_moments

         use parameters_numerical, only: autostop

         implicit none

         !-------------------------------------------------------------------------

         ! Track code 
         if (debug) write (*, *) 'read_diagnostics_parameters::read_input_file'


         call read_namelist_diagnostics (nwrite, navg, nsave, nc_mult, &
                        save_for_restart, write_all, write_all_time_traces, &
                        write_all_spectra_kxkyz, write_all_spectra_kxky, &
                        write_all_velocity_space, write_all_potential, &
                        write_all_omega, write_all_distribution, write_all_fluxes, &
                        write_all_moments)

         call read_namelist_diagnostics_potential (write_all_potential, &
                        write_all_time_traces, write_all_spectra_kxkyz, &
                        write_all_spectra_kxky, &
                        write_all_potential_time_traces, write_all_potential_spectra, &
                        write_phi2_vs_time, write_apar2_vs_time, &
                        write_bpar2_vs_time, write_phi_vs_kxkyz, &
                        write_apar_vs_kxkyz, write_bpar_vs_kxkyz, &
                        write_phi2_vs_kxky, write_apar2_vs_kxky, & 
                        write_bpar2_vs_kxky)
         
         call read_namelist_diagnostics_omega (write_all_omega, &
                        write_omega_vs_kxky, write_omega_avg_vs_kxky)

         call read_namelist_diagnostics_distribution (write_all_distribution, &
                        write_all_spectra_kxkyz, &
                        write_all_velocity_space, &
                        write_g2_vs_vpamus, write_g2_vs_zvpas, &
                        write_g2_vs_zmus, write_g2_vs_kxkyzs, &
                        write_g2_vs_zvpamus, &
                        write_distribution_g, write_distribution_h, &
                        write_distribution_f)
         
         call read_namelist_diagnostics_fluxes (write_all_fluxes, &
                        write_all_time_traces, write_all_spectra_kxkyz, &
                        write_all_spectra_kxky, &
                        flux_norm, write_fluxes_vs_time, &
                        write_radial_fluxes, write_fluxes_kxkyz, &
                        write_fluxes_kxky)

         call read_namelist_diagnostics_moments (write_all_moments, &
                        write_moments, write_radial_moments)
            
         !-------------------------------------------------------------------
         
         ! If <save_for_restart> = False then we need <nsave> = -1
         if (.not. save_for_restart) nsave = -1

         ! For nonlinear simulations, don't stop automatically 
         ! and if we want to <autostop> we need the <omega_vs_tkykx> 
         ! which contains <omega_vs_kykx> over the last <navg> time steps
         if (autostop) write_omega_avg_vs_kxky = .true.

      end subroutine read_input_file 
   
      !**********************************************************************
      !                         BROADCAST PARAMETERS                        !
      !**********************************************************************
      subroutine broadcast_parameters
      
         use mp, only: broadcast
      
         implicit none

         !-------------------------------------------------------------------------
         
         ! Track code 
         if (debug) write (*, *) 'read_diagnostics_parameters::broadcast_parameters'
       
         ! Writing options 
         call broadcast(nwrite)
         call broadcast(navg)
         call broadcast(nsave)
         call broadcast(nc_mult)
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
         call broadcast(write_phi_vs_kxkyz)
         call broadcast(write_apar_vs_kxkyz)
         call broadcast(write_bpar_vs_kxkyz)
         call broadcast(write_phi2_vs_kxky)
         call broadcast(write_apar2_vs_kxky)
         call broadcast(write_bpar2_vs_kxky)
         call broadcast(write_moments)
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

   end subroutine read_parameters_diagnostics

end module parameters_diagnostics
