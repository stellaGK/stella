# include "define.inc"
!###############################################################################
!                   SAVE DISTRIBUTION FUNCTION FOR A RESTART                    
!###############################################################################
! This module saves the distribution function g(mu,vpa,ikxkyzs) to a netcdf
! file, as well as the dimensions (kx,ky,z,nu,vpa,s,tube). This allows us to
! restart the simulation in the future, to continue the time evolution.
!###############################################################################
module save_stella_for_restart

   ! Import mpi
   use mp, only: mp_comm
   use mp, only: mp_info

   ! Import netcdf modules
#ifdef NETCDF
   use netcdf, only: NF90_NOWRITE, NF90_CLOBBER, NF90_NOERR
   use netcdf, only: nf90_create, nf90_open, nf90_sync, nf90_close
   use netcdf, only: nf90_def_dim, nf90_def_var, nf90_enddef
   use netcdf, only: nf90_put_var, nf90_get_var, nf90_strerror
   use netcdf, only: nf90_inq_dimid, nf90_inquire_dimension
   use netcdf, only: nf90_inq_varid, nf90_inquire_variable
   use netcdf, only: nf90_int
   use netcdf_utils, only: get_netcdf_code_precision
   use netcdf_utils, only: check_netcdf_file_precision
   use netcdf_utils, only: netcdf_error
   use netcdf_utils, only: netcdf_real, kind_nf
#endif

   implicit none

   ! Make the routines available to other modules
   public :: stella_restore
   public :: save_stella_data_for_restart
   public :: save_many
   public :: init_dt
   public :: init_tstart
   public :: init_save
 
   ! Note that <save_many> can be set in the input file in <restart_options>
   ! but it is only read if <save_for_restart> = True. If <save_many> = True,
   ! the distribution function is saved to a netcdf file for each processor.
   ! If instead <save_many> = .false., it is saved to a single netcdf file.
   logical :: save_many = .true.

   private
   
   ! Name of the restart file
   character(300), save :: restart_file

   ! Initialise local arrays
#ifdef NETCDF
   integer(kind_nf) :: ncid, zedid, vpaid, gloid, gvmuloid, kyid, kxid, muid, tubeid
   integer(kind_nf) :: krookr_id, krooki_id, projr_id, proji_id
   integer(kind_nf) :: phiprojr_id, phiproji_id
   integer(kind_nf) :: t0id, gr_id, gi_id, delt0id, istep0id
   integer(kind_nf) :: intkrook_id, intproj_id;
   integer(kind_nf) :: shift_id
   logical :: initialised_restart_module = .false.
#endif
   logical :: intialised_parallel_netcdf = .false.
   logical :: parallel_netcdf

contains


!###############################################################################
!############################### SAVE DISTRIBUTION #############################
!###############################################################################

   !****************************************************************************
   !                 Save distribution function to a netcdf file                
   !****************************************************************************
   subroutine save_stella_data_for_restart(istep0, t0, delt0, istatus, exit_in)

#ifdef NETCDF_PARALLEL
      use save_stella_for_restart_parallel_netcdf, only: save_stella_data_for_restart_to_a_single_file
#endif
      
      implicit none

      ! Arguments
      real, intent(in) :: t0, delt0
      integer, intent(in) :: istep0
      integer, intent(out) :: istatus
      logical, intent(in), optional :: exit_in

      !-------------------------------------------------------------------------

      ! Initialise the restart module
      if (.not. intialised_parallel_netcdf) then

         ! Only initialise once
         intialised_parallel_netcdf = .true.
         
         ! The distribution function can be saved to a single restart file, or it
         ! can be saved to multiple restart file, using one file for each processor.
         ! To save to a single restart file, parallel netcdf needs to be loaded.
         ! The save to multiple restart files, set <save_many> = True.
         parallel_netcdf = .false.
#ifdef NETCDF_PARALLEL
         parallel_netcdf = .true.
#endif
         
         ! The <save_many> variable is an input variable that can be set by the user
         ! However, if parallel netcdf is not loaded, we need to save to multiple netcdf files
         if (.not. parallel_netcdf) save_many = .true.
         
      end if
      
      ! Save to multiple netcdf files or a single netcdf file
#ifdef NETCDF
      if (save_many) call save_stella_data_for_restart_to_multiple_files(istep0, t0, delt0, istatus, exit_in)
#ifdef NETCDF_PARALLEL
      if (.not. save_many) call save_stella_data_for_restart_to_a_single_file(istep0, t0, delt0, istatus, exit_in, restart_file)
#endif
#else
      write(*,*) 'Could not save restart data to a netcdf file, since netcdf is not loaded.'
#endif
      
   end subroutine save_stella_data_for_restart

#ifdef NETCDF
   !****************************************************************************
   !    Save distribution function to a netcdf file (one file per processor)    
   !****************************************************************************
   subroutine save_stella_data_for_restart_to_multiple_files(istep0, t0, delt0, istatus, exit_in)
      
      ! Must include kxkyz_layout_type here to avoid obscure bomb while compiling
      ! diagnostics.f90 (which uses this module) with the Compaq F90 compiler:
      use parallelisation_layouts, only: kxkyz_lo, vmu_lo
      use common_types, only: kxkyz_layout_type
      use mp, only: iproc
      
      ! Grids
      use grids_z, only: nzgrid, ntubes
      use grids_velocity, only: nvpa, nmu
      use grids_kxky, only: naky, nakx
      
      ! Flow shear
      use arrays, only: shift_state
      
      ! Make sure we have the distribution function vs (nvpa, nmu, -kxkyz-layout-)
      use redistribute, only: scatter
      use initialise_redistribute, only: kxkyz2vmu
      use arrays_distribution_function, only: gnew
      use arrays_distribution_function, only: gvmu

      implicit none

      ! Arguments
      real, intent(in) :: t0, delt0
      integer, intent(in) :: istep0
      integer, intent(out) :: istatus
      logical, intent(in), optional :: exit_in
      
      ! Local variables
      character(306) :: path_netcdf_file_per_proc
      character(10) :: suffix
      integer :: n_elements, nvmulo_elements
      logical :: has_vmulo
      logical :: exit

      ! Temporary arrays
      real, allocatable, dimension(:, :, :) :: g_temp
      
      !-------------------------------------------------------------------------
      
      ! Error status and exit status
      istatus = 0
      if (present(exit_in)) then
         exit = exit_in
      else
         exit = .false.
      end if
      
      ! Make sure we have the distribution function vs (nvpa, nmu, -kxkyz-layout-)
      call scatter(kxkyz2vmu, gnew, gvmu)

      ! Parallelisation
      n_elements = kxkyz_lo%ulim_proc - kxkyz_lo%llim_proc + 1
      nvmulo_elements = vmu_lo%ulim_proc - vmu_lo%llim_proc + 1
      if (n_elements <= 0) return

      ! Create netcdf file for each vmulo element
      has_vmulo = nvmulo_elements > 0 .or. .not. save_many

      !-------------------------------------------------------------------------
      !                       Initialise the netcdf file                        
      !-------------------------------------------------------------------------
      if (.not. initialised_restart_module) then

         ! Only initialise once
         initialised_restart_module = .true.

         ! Path of the netcdf file for each processor
         write (suffix, '(a1,i0)') '.', iproc
         path_netcdf_file_per_proc = trim(restart_file)
         path_netcdf_file_per_proc = trim(trim(path_netcdf_file_per_proc)//adjustl(suffix))
         
         ! Create a netcdf file for each processor
         istatus = nf90_create(path_netcdf_file_per_proc, NF90_CLOBBER, ncid)
         if (istatus /= NF90_NOERR) call process_nf90_error("nf90_create error: ", istatus)

         ! Save the dimensions
         if (n_elements > 0) then
         
            ! Save the (kx,ky,z,tube,nu,vpa) dimensions
            istatus = nf90_def_dim(ncid, "akx", nakx, kxid)
            if (istatus /= NF90_NOERR) call process_nf90_error("nf90_def_dim akx error: ", istatus)
            istatus = nf90_def_dim(ncid, "aky", naky, kyid)
            if (istatus /= NF90_NOERR) call process_nf90_error("nf90_def_dim aky error: ", istatus)
            istatus = nf90_def_dim(ncid, "zed", 2 * nzgrid + 1, zedid)
            if (istatus /= NF90_NOERR) call process_nf90_error("nf90_def_dim zed error: ", istatus)
            istatus = nf90_def_dim(ncid, "tube", ntubes, tubeid)
            if (istatus /= NF90_NOERR) call process_nf90_error("nf90_def_dim tube error: ", istatus)
            istatus = nf90_def_dim(ncid, "vpa", nvpa, vpaid)
            if (istatus /= NF90_NOERR) call process_nf90_error("nf90_def_dim vpa error: ", istatus)
            istatus = nf90_def_dim(ncid, "mu", nmu, muid)
            if (istatus /= NF90_NOERR) call process_nf90_error("nf90_def_dim mu error: ", istatus)
            
            ! Save the number of ikxkyzs points
            istatus = nf90_def_dim(ncid, "glo", n_elements, gloid)
            if (istatus /= NF90_NOERR) call process_nf90_error("nf90_def_dim glo error: ", istatus)

            ! Save the number of ivpamus points
            if (nvmulo_elements > 0) then
               istatus = nf90_def_dim(ncid, "gvmulo", nvmulo_elements, gvmuloid)
               if (istatus /= NF90_NOERR) call process_nf90_error("nf90_def_dim gvmulo error: ", istatus)
            end if

         end if

         ! At initialisation there the variables swill set to zero
         ! The zero is chosen to be code precision instead of absolute zero
         if (netcdf_real == 0) netcdf_real = get_netcdf_code_precision()

         ! Save the time, istep and delt variables
         istatus = nf90_def_var(ncid, "t0", netcdf_real, t0id)
         if (istatus /= NF90_NOERR) call process_nf90_error("nf90_def_var t0 error: ", istatus)
         istatus = nf90_def_var(ncid, "istep0", nf90_int, istep0id)
         if (istatus /= NF90_NOERR) call process_nf90_error("nf90_def_var istep0 error: ", istatus)
         istatus = nf90_def_var(ncid, "delt0", netcdf_real, delt0id)
         if (istatus /= NF90_NOERR) call process_nf90_error("nf90_def_var delt0 error: ", istatus)

         ! Save the variables
         if (n_elements > 0) then
         
            ! Save the real and imaginary part of the distribution function vs (vpa, mu, ikxkyzs)
            istatus = nf90_def_var(ncid, "gr", netcdf_real, (/vpaid, muid, gloid/), gr_id)
            if (istatus /= NF90_NOERR) call process_nf90_error("nf90_def_var g error: ", istatus)
            istatus = nf90_def_var(ncid, "gi", netcdf_real, (/vpaid, muid, gloid/), gi_id)
            if (istatus /= NF90_NOERR) call process_nf90_error("nf90_def_var g error: ", istatus)

            ! Flow shear
            istatus = nf90_def_var(ncid, "shiftstate", netcdf_real, (/kyid/), shift_id)
            if (istatus /= NF90_NOERR) call process_nf90_error("nf90_def_var shiftstate error: ", istatus)
            
            ! Radial variation and Sources
            call save_radial_variation_to_netcdf_def_var()

         end if

         ! Finished defining all dimensions and variables in the netcdf file
         istatus = nf90_enddef(ncid)
         if (istatus /= NF90_NOERR) call process_nf90_error("nf90_enddef error: ", istatus)

      end if

      !-------------------------------------------------------------------------
      !                          Write the netcdf file                          
      !-------------------------------------------------------------------------

      ! Write the actual time, istep and delt variables to the netcdf file
      istatus = nf90_put_var(ncid, t0id, t0)
      if (istatus /= NF90_NOERR) call process_nf90_error("nf90_put_var t0 error: ", istatus)
      istatus = nf90_put_var(ncid, istep0id, istep0)
      if (istatus /= NF90_NOERR) call process_nf90_error("nf90_put_var istep error: ", istatus)
      istatus = nf90_put_var(ncid, delt0id, delt0)
      if (istatus /= NF90_NOERR) call process_nf90_error("nf90_put_var delt0 error: ", istatus)

      ! Save the distribution function to the netcdf file
      if (n_elements > 0) then

         ! Use a temporary array for the real and imaginary parts of g(vpa,mu,ikxkyzs)
         if (.not. allocated(g_temp)) allocate (g_temp(nvpa, nmu, kxkyz_lo%llim_proc:kxkyz_lo%ulim_alloc))

         ! Save the real part of the distribution function
         g_temp = real(gvmu)
         istatus = nf90_put_var(ncid, gr_id, g_temp)
         if (istatus /= NF90_NOERR) call netcdf_error(istatus, ncid, gr_id)

         ! Save the imaginary part of the distribution function
         g_temp = aimag(gvmu)
         istatus = nf90_put_var(ncid, gi_id, g_temp)
         if (istatus /= NF90_NOERR) call netcdf_error(istatus, ncid, gi_id)
         
         ! Deallocate the temporary array
         if (allocated(g_temp)) deallocate (g_temp)

         ! Flow shear
         istatus = nf90_put_var(ncid, shift_id, shift_state)
         if (istatus /= NF90_NOERR) call netcdf_error(istatus, ncid, shift_id)
         
         ! Radial variation and sources
         call save_radial_variation_to_netcdf_put_var()
         
      end if

      ! If we make a clean exit of stella, then close the netcdf file
      if (exit) then
         istatus = nf90_close(ncid)
         if (istatus /= NF90_NOERR) call netcdf_error(istatus, message='nf90_close error')
         
      ! Otherwise, sync the netcdf file
      else
         istatus = nf90_sync(ncid)
         if (istatus /= NF90_NOERR) call netcdf_error(istatus, message='nf90_sync error')
      end if
      
   contains
   
      !------------- Radial variation and Sources: Define variables ------------
      subroutine save_radial_variation_to_netcdf_def_var
      
         ! Radial variation and sources
         use gk_sources, only: include_qn_source
         use gk_sources, only: source_option_krook
         use gk_sources, only: source_option_projection
         use gk_sources, only: source_option_switch
         
         implicit none

         !----------------------------------------------------------------------
      
         ! Radial variation and sources
         if (source_option_switch == source_option_krook .and. has_vmulo) then
            istatus = nf90_def_var(ncid, "intkrook", netcdf_real, intkrook_id)
            if (istatus /= NF90_NOERR) call process_nf90_error("nf90_def_var intkrook error: ", istatus)
            istatus = nf90_def_var(ncid, "krookr", netcdf_real, (/kxid, zedid, tubeid, gvmuloid/), krookr_id)
            if (istatus /= NF90_NOERR) call process_nf90_error("nf90_def_var krookr error: ", istatus)
            istatus = nf90_def_var(ncid, "krooki", netcdf_real, (/kxid, zedid, tubeid, gvmuloid/), krooki_id)
            if (istatus /= NF90_NOERR) call process_nf90_error("nf90_def_var krooki error: ", istatus)
         end if

         ! Radial variation and sources
         if (include_qn_source .and. iproc == 0) then
            istatus = nf90_def_var(ncid, "phiprojr", netcdf_real, (/kxid, zedid, tubeid/), phiprojr_id)
            if (istatus /= NF90_NOERR) call process_nf90_error("nf90_def_var phiprojr error: ", istatus)
            istatus = nf90_def_var(ncid, "phiproji", netcdf_real, (/kxid, zedid, tubeid/), phiproji_id)
            if (istatus /= NF90_NOERR) call process_nf90_error("nf90_def_var phiproji error: ", istatus)
         end if

         ! Radial variation and sources
         if (source_option_switch == source_option_projection .and. has_vmulo) then
            istatus = nf90_def_var(ncid, "intproj", netcdf_real, intproj_id)
            if (istatus /= NF90_NOERR) call process_nf90_error("nf90_def_var intproj error: ", istatus)
            istatus = nf90_def_var(ncid, "projr", netcdf_real, (/kxid, zedid, tubeid, gvmuloid/), projr_id)
            if (istatus /= NF90_NOERR) call process_nf90_error("nf90_def_var projr error: ", istatus)
            istatus = nf90_def_var(ncid, "proji", netcdf_real, (/kxid, zedid, tubeid, gvmuloid/), proji_id)
            if (istatus /= NF90_NOERR) call process_nf90_error("nf90_def_var proji error: ", istatus)
         end if
      
      end subroutine save_radial_variation_to_netcdf_def_var
      
      !-------------- Radial variation and Sources: Put variables --------------
      subroutine save_radial_variation_to_netcdf_put_var
      
         ! Radial variation and sources
         use gk_sources, only: include_qn_source
         use gk_sources, only: source_option_krook, source_option_projection
         use gk_sources, only: source_option_switch, int_krook, int_proj
         use arrays_distribution_function, only: g_krook, g_proj
         use arrays_fields, only: phi_proj
         
         implicit none
         
         real, allocatable, dimension(:, :, :, :) :: krook_temp
         real, allocatable, dimension(:, :, :, :) :: proj_temp
         real, allocatable, dimension(:, :, :) :: phiproj_temp

         !----------------------------------------------------------------------

         ! Radial variation and sources
         if (source_option_switch == source_option_krook .and. has_vmulo) then
            istatus = nf90_put_var(ncid, intkrook_id, int_krook)
            if (istatus /= NF90_NOERR) call process_nf90_error("nf90_put_var int_krook error: ", istatus)
            if (.not. allocated(krook_temp)) allocate (krook_temp(nakx, -nzgrid:nzgrid, ntubes, vmu_lo%llim_proc:vmu_lo%ulim_alloc))
            krook_temp = real(g_krook)
            istatus = nf90_put_var(ncid, krookr_id, krook_temp)
            if (istatus /= NF90_NOERR) call netcdf_error(istatus, ncid, krookr_id)
            krook_temp = aimag(g_krook)
            istatus = nf90_put_var(ncid, krooki_id, krook_temp)
            if (istatus /= NF90_NOERR) call netcdf_error(istatus, ncid, krooki_id)
            if (allocated(krook_temp)) deallocate (krook_temp)
         end if

         ! Radial variation and sources
         if (source_option_switch == source_option_projection .and. has_vmulo) then
            istatus = nf90_put_var(ncid, intproj_id, int_proj)
            if (istatus /= NF90_NOERR) call process_nf90_error("nf90_put_var int_proj error: ", istatus)
            if (.not. allocated(proj_temp)) allocate (proj_temp(nakx, -nzgrid:nzgrid, ntubes, vmu_lo%llim_proc:vmu_lo%ulim_alloc))
            proj_temp = real(g_proj)
            istatus = nf90_put_var(ncid, projr_id, proj_temp)
            if (istatus /= NF90_NOERR) call netcdf_error(istatus, ncid, projr_id)
            proj_temp = aimag(g_proj)
            istatus = nf90_put_var(ncid, proji_id, proj_temp)
            if (istatus /= NF90_NOERR) call netcdf_error(istatus, ncid, proji_id)
            if (allocated(proj_temp)) deallocate (proj_temp)
         end if

         ! Radial variation and sources
         if (include_qn_source .and. iproc == 0) then
            if (.not. allocated(phiproj_temp)) allocate (phiproj_temp(nakx, -nzgrid:nzgrid, ntubes))
            phiproj_temp = real(phi_proj)
            istatus = nf90_put_var(ncid, phiprojr_id, phiproj_temp)
            if (istatus /= NF90_NOERR) call netcdf_error(istatus, ncid, phiprojr_id)
            phiproj_temp = aimag(phi_proj)
            istatus = nf90_put_var(ncid, phiproji_id, phiproj_temp)
            if (istatus /= NF90_NOERR) call netcdf_error(istatus, ncid, phiproji_id)
            if (allocated(phiproj_temp)) deallocate (phiproj_temp)
         end if
      
      end subroutine save_radial_variation_to_netcdf_put_var

   end subroutine save_stella_data_for_restart_to_multiple_files
#endif

!###############################################################################
!############################### READ DISTRIBUTION #############################
!###############################################################################

   !****************************************************************************
   !                Read distribution function from a netcdf file               
   !****************************************************************************
   subroutine stella_restore(g, scale, istatus)
   
      use parallelisation_layouts, only: kxkyz_lo
#ifdef NETCDF_PARALLEL
      use save_stella_for_restart_parallel_netcdf, only: read_stella_data_for_restart_from_single_file
#endif

      implicit none

      ! Arguments
      complex, dimension(:, :, kxkyz_lo%llim_proc:), intent(in out) :: g
      real, intent(in) :: scale
      integer, intent(out) :: istatus

      !-------------------------------------------------------------------------
      
      ! Initialise the restart module
      if (.not. intialised_parallel_netcdf) then

         ! Only initialise once
         intialised_parallel_netcdf = .true.
         
         ! The distribution function can be saved to a single restart file, or it
         ! can be saved to multiple restart file, using one file for each processor.
         ! To save to a single restart file, parallel netcdf needs to be loaded.
         ! The save to multiple restart files, set <save_many> = True.
         parallel_netcdf = .false.
#ifdef NETCDF_PARALLEL
         parallel_netcdf = .true.
#endif
         
         ! The <save_many> variable is an input variable that can be set by the user
         ! However, if parallel netcdf is not loaded, we need to save to multiple netcdf files
         if (.not. parallel_netcdf) save_many = .true.
         
      end if
      
      ! Read from multiple netcdf files or a single netcdf file
#ifdef NETCDF
      if (save_many) call read_stella_data_for_restart_from_multiple_files(g, scale, istatus)
#ifdef NETCDF_PARALLEL
      if (.not. save_many) call read_stella_data_for_restart_from_single_file(g, scale, istatus, restart_file)
#endif
#else
      write(*,*) 'Could not read restart data from a netcdf file, since netcdf is not loaded.'
#endif
      
   end subroutine stella_restore

#ifdef NETCDF

   !****************************************************************************
   !            Read distribution function from multiple netcdf files           
   !****************************************************************************
   subroutine read_stella_data_for_restart_from_multiple_files(g, scale, istatus)
   
      ! Parallelisation
      use mp, only: iproc, broadcast
      use parallelisation_layouts, only: kxkyz_lo, vmu_lo
      
      ! Grids
      use grids_kxky, only: naky, nakx
      use grids_z, only: nzgrid, ntubes
      use grids_velocity, only: nvpa, nmu
      
      ! Flow shear
      use arrays, only: shift_state

      implicit none

      ! Arguments
      complex, dimension(:, :, kxkyz_lo%llim_proc:), intent(out) :: g
      real, intent(in) :: scale
      integer, intent(out) :: istatus
      
      ! Local variables
      character(306) :: path_netcdf_file_per_proc
      character(10) :: suffix
      integer :: i, n_elements, nvmulo_elements
      logical :: has_vmulo

      ! Temporary arrays
      real, allocatable, dimension(:, :, :) :: g_real_temp, g_imag_temp

      !-------------------------------------------------------------------------
       
      ! Number of ikxkyzs points and imuvpas points
      n_elements = kxkyz_lo%ulim_proc - kxkyz_lo%llim_proc + 1
      nvmulo_elements = vmu_lo%ulim_proc - vmu_lo%llim_proc + 1

      ! Don't read if something went wrong with the parallelisation
      if (n_elements <= 0) return

      ! Check if we parallelised over velocity space
      has_vmulo = nvmulo_elements > 0

      !-------------------------------------------------------------------------
      !                       Initialise the netcdf file                        
      !-------------------------------------------------------------------------
      if (.not. initialised_restart_module) then
      
         ! Path of the netcdf file for each processor
         write (suffix, '(a1,i0)') '.', iproc
         path_netcdf_file_per_proc = trim(restart_file) 
         path_netcdf_file_per_proc = trim(trim(path_netcdf_file_per_proc)//adjustl(suffix))
         
         ! Open the netcdf file for each processor
         istatus = nf90_open(path_netcdf_file_per_proc, NF90_NOWRITE, ncid)
         if (istatus /= NF90_NOERR) call netcdf_error(istatus, file=path_netcdf_file_per_proc, abort=.true.)

         ! Check code precision
         if (netcdf_real == 0) netcdf_real = get_netcdf_code_precision()
         call check_netcdf_file_precision(ncid)

         ! Inquire about the netcdf dimensions id
         istatus = nf90_inq_dimid(ncid, "akx", kxid)
         if (istatus /= NF90_NOERR) call netcdf_error(istatus, dim='akx')
         istatus = nf90_inq_dimid(ncid, "aky", kyid)
         if (istatus /= NF90_NOERR) call netcdf_error(istatus, dim='aky')
         istatus = nf90_inq_dimid(ncid, "zed", zedid)
         if (istatus /= NF90_NOERR) call netcdf_error(istatus, dim='zed')
         istatus = nf90_inq_dimid(ncid, "tube", tubeid)
         if (istatus /= NF90_NOERR) call netcdf_error(istatus, dim='tube')
         
         ! Number of ikxkyzs points
         istatus = nf90_inq_dimid(ncid, "glo", gloid)
         if (istatus /= NF90_NOERR) call netcdf_error(istatus, dim='glo')
         
         ! Nmber of ivpamus points
         if (has_vmulo) then
            istatus = nf90_inq_dimid(ncid, "gvmulo", gvmuloid)
            if (istatus /= NF90_NOERR) call netcdf_error(istatus, dim='gvmulo')
         end if

         ! Check that the dimensions in the netcdf file match those specified in the input file
         istatus = nf90_inquire_dimension(ncid, kxid, len=i)
         if (istatus /= NF90_NOERR) call netcdf_error(istatus, ncid, dimid=kxid)
         if (i /= nakx) write (*, *) 'Restart error: nakx=? ', i, ' : ', nakx, ' : ', iproc
         istatus = nf90_inquire_dimension(ncid, kyid, len=i)
         if (istatus /= NF90_NOERR) call netcdf_error(istatus, ncid, dimid=kyid)
         if (i /= naky) write (*, *) 'Restart error: naky=? ', i, ' : ', naky, ' : ', iproc
         istatus = nf90_inquire_dimension(ncid, zedid, len=i)
         if (istatus /= NF90_NOERR) call netcdf_error(istatus, ncid, dimid=zedid)
         if (i /= 2 * nzgrid + 1) write (*, *) 'Restart error: nzgrid=? ', i, ' : ', nzgrid, ' : ', iproc
         istatus = nf90_inquire_dimension(ncid, tubeid, len=i)
         if (istatus /= NF90_NOERR) call netcdf_error(istatus, ncid, dimid=tubeid)
         if (i /= ntubes) write (*, *) 'Restart error: ntubes=? ', i, ' : ', ntubes, ' : ', iproc
         istatus = nf90_inquire_dimension(ncid, gloid, len=i)
         if (istatus /= NF90_NOERR) call netcdf_error(istatus, ncid, dimid=gloid)
         if (i /= n_elements) write (*, *) 'Restart error: glo=? ', i, ' : ', iproc
         if (has_vmulo) then
            istatus = nf90_inquire_dimension(ncid, gvmuloid, len=i)
            if (istatus /= NF90_NOERR) call netcdf_error(istatus, ncid, dimid=gvmuloid)
            if (i /= nvmulo_elements) write (*, *) 'Restart error: gvmulo=? ', i, ' : ', iproc
         end if

         ! Distribution function ids
         istatus = nf90_inq_varid(ncid, "gr", gr_id)
         if (istatus /= NF90_NOERR) call netcdf_error(istatus, var='gr')
         istatus = nf90_inq_varid(ncid, "gi", gi_id)
         if (istatus /= NF90_NOERR) call netcdf_error(istatus, var='gi')
         
         ! Flow shear ids
         istatus = nf90_inq_varid(ncid, "shiftstate", shift_id)
         if (istatus /= NF90_NOERR) call netcdf_error(istatus, var='shiftstate')

         ! Radial variation and Sources ids
         call read_radial_variation_from_netcdf_variable_ids()
   
      end if

      !-------------------------------------------------------------------------
      !                          Read the netcdf file                          
      !-------------------------------------------------------------------------
      
      ! Allocate the arrays for the distribution function and set them to zero
      if (.not. allocated(g_real_temp)) allocate (g_real_temp(nvpa, nmu, kxkyz_lo%llim_proc:kxkyz_lo%ulim_alloc))
      if (.not. allocated(g_imag_temp)) allocate (g_imag_temp(nvpa, nmu, kxkyz_lo%llim_proc:kxkyz_lo%ulim_alloc))
      g_real_temp = 0.; g_imag_temp = 0.
      
      ! Read the real part of the distribution function
      istatus = nf90_get_var(ncid, gr_id, g_real_temp)
      if (istatus /= NF90_NOERR) call netcdf_error(istatus, ncid, gr_id)

      ! Read the imaginary part of the distribution function
      istatus = nf90_get_var(ncid, gi_id, g_imag_temp)
      if (istatus /= NF90_NOERR) call netcdf_error(istatus, ncid, gi_id)

      ! Construct the distribution function vs (vpa, mu, ikxkyzs)
      g = cmplx(g_real_temp, g_imag_temp)
      
      ! Deallocate temporary arrays
      if (allocated(g_real_temp)) deallocate (g_real_temp)
      if (allocated(g_imag_temp)) deallocate (g_imag_temp)

      ! Flow shear
      if (.not. allocated(shift_state)) allocate (shift_state(naky))
      istatus = nf90_get_var(ncid, shift_id, shift_state)
      if (istatus /= NF90_NOERR) call netcdf_error(istatus, ncid, shift_id)
      
      ! Radial variation and sources
      call read_radial_variation_from_netcdf_variables(scale)

      ! When restarting linear simulations, we can rescale the potential
      ! to avoid running into infinity numbers
      if (scale > 0.) g = g * scale
      
   contains
   
      !-------------- Radial variation and Sources: Variable ids --------------
      subroutine read_radial_variation_from_netcdf_variable_ids
      
         ! Radial variation and sources
         use gk_sources, only: include_qn_source
         use gk_sources, only: source_option_krook
         use gk_sources, only: source_option_projection
         use gk_sources, only: source_option_switch
         
         implicit none

         !----------------------------------------------------------------------
      
         if (source_option_switch == source_option_krook .and. has_vmulo) then
            istatus = nf90_inq_varid(ncid, "intkrook", intkrook_id)
            if (istatus /= NF90_NOERR) call netcdf_error(istatus, var='intkrook')
            istatus = nf90_inq_varid(ncid, "krookr", krookr_id)
            if (istatus /= NF90_NOERR) call netcdf_error(istatus, var='krookr')
            istatus = nf90_inq_varid(ncid, "krooki", krooki_id)
            if (istatus /= NF90_NOERR) call netcdf_error(istatus, var='krooki')
         end if

         ! Radial variation and Sources ids
         if (source_option_switch == source_option_projection .and. has_vmulo) then
            istatus = nf90_inq_varid(ncid, "intproj", intproj_id)
            if (istatus /= NF90_NOERR) call netcdf_error(istatus, var='intproj')
            istatus = nf90_inq_varid(ncid, "projr", projr_id)
            if (istatus /= NF90_NOERR) call netcdf_error(istatus, var='projr')
            istatus = nf90_inq_varid(ncid, "proji", proji_id)
            if (istatus /= NF90_NOERR) call netcdf_error(istatus, var='proji')
         end if

         ! Radial variation and Sources ids
         if (include_qn_source .and. iproc == 0) then
            istatus = nf90_inq_varid(ncid, "phiprojr", phiprojr_id)
            if (istatus /= NF90_NOERR) call netcdf_error(istatus, var='phiprojr')
            istatus = nf90_inq_varid(ncid, "phiproji", phiproji_id)
            if (istatus /= NF90_NOERR) call netcdf_error(istatus, var='phiproji')
         end if
      
      end subroutine read_radial_variation_from_netcdf_variable_ids
      
      !-------------- Radial variation and Sources: Read variables -------------
      subroutine read_radial_variation_from_netcdf_variables(scale)
      
         ! Radial variation and sources
         use gk_sources, only: include_qn_source
         use gk_sources, only: source_option_krook, source_option_projection
         use gk_sources, only: source_option_switch, int_krook, int_proj
         use arrays_distribution_function, only: g_krook, g_proj
         use arrays_fields, only: phi_proj
         
         implicit none
         
         ! Arguments
         real, intent(in) :: scale
         
         ! Temporary arrays
         real, allocatable, dimension(:, :, :, :) :: ktmpr, ktmpi
         real, allocatable, dimension(:, :, :, :) :: ptmpr, ptmpi
         real, allocatable, dimension(:, :, :) :: pptmpr, pptmpi

         !----------------------------------------------------------------------

         ! Radial variation and sources
         if (source_option_switch == source_option_krook .and. has_vmulo) then
            istatus = nf90_get_var(ncid, intkrook_id, int_krook)
            if (istatus /= NF90_NOERR) call netcdf_error(istatus, ncid, intkrook_id)
            if (.not. allocated(ktmpr)) allocate (ktmpr(nakx, -nzgrid:nzgrid, ntubes, vmu_lo%llim_proc:vmu_lo%ulim_alloc))
            if (.not. allocated(ktmpi)) allocate (ktmpi(nakx, -nzgrid:nzgrid, ntubes, vmu_lo%llim_proc:vmu_lo%ulim_alloc))
            ktmpr = 0.; ktmpi = 0.
            istatus = nf90_get_var(ncid, krookr_id, ktmpr)
            if (istatus /= NF90_NOERR) call netcdf_error(istatus, ncid, krookr_id)
            istatus = nf90_get_var(ncid, krooki_id, ktmpi)
            if (istatus /= NF90_NOERR) call netcdf_error(istatus, ncid, krooki_id)
            g_krook = cmplx(ktmpr, ktmpi)
            if (allocated(ktmpr)) deallocate (ktmpr)
            if (allocated(ktmpi)) deallocate (ktmpi)
         end if

         ! Radial variation and sources
         if (source_option_switch == source_option_projection .and. has_vmulo) then
            istatus = nf90_get_var(ncid, intproj_id, int_proj)
            if (istatus /= NF90_NOERR) call netcdf_error(istatus, ncid, intproj_id)
            if (.not. allocated(ptmpr)) allocate (ptmpr(nakx, -nzgrid:nzgrid, ntubes, vmu_lo%llim_proc:vmu_lo%ulim_alloc))
            if (.not. allocated(ptmpi)) allocate (ptmpi(nakx, -nzgrid:nzgrid, ntubes, vmu_lo%llim_proc:vmu_lo%ulim_alloc))
            ptmpr = 0.; ptmpi = 0.
            istatus = nf90_get_var(ncid, projr_id, ptmpr)
            if (istatus /= NF90_NOERR) call netcdf_error(istatus, ncid, projr_id)
            istatus = nf90_get_var(ncid, proji_id, ptmpi)
            if (istatus /= NF90_NOERR) call netcdf_error(istatus, ncid, proji_id)
            g_proj = cmplx(ptmpr, ptmpi)
            if (allocated(ptmpr)) deallocate (ptmpr)
            if (allocated(ptmpi)) deallocate (ptmpi)
         end if

         ! Radial variation and sources
         if (include_qn_source .and. iproc == 0) then
            if (.not. allocated(pptmpr)) allocate (pptmpr(nakx, -nzgrid:nzgrid, ntubes))
            if (.not. allocated(pptmpi)) allocate (pptmpi(nakx, -nzgrid:nzgrid, ntubes))
            pptmpr = 0.; pptmpi = 0.
            istatus = nf90_get_var(ncid, phiprojr_id, pptmpr)
            if (istatus /= NF90_NOERR) call netcdf_error(istatus, ncid, phiprojr_id)
            istatus = nf90_get_var(ncid, phiproji_id, pptmpi)
            if (istatus /= NF90_NOERR) call netcdf_error(istatus, ncid, phiproji_id)
            phi_proj = cmplx(pptmpr, pptmpi)
            if (allocated(pptmpr)) deallocate (pptmpr)
            if (allocated(pptmpi)) deallocate (pptmpi)
         end if
         
         ! When restarting linear simulations, we can rescale the potential
         if (scale > 0.) then
            if (source_option_switch == source_option_krook) g_krook = g_krook * scale
            if (source_option_switch == source_option_projection) g_proj = g_proj * scale
         end if
            
         ! Broadcast phi_proj
         if (include_qn_source) call broadcast(phi_proj)
            
      end subroutine read_radial_variation_from_netcdf_variables

   end subroutine read_stella_data_for_restart_from_multiple_files
#endif
   
   !****************************************************************************
   !                                      Title
   !****************************************************************************
   subroutine init_save(file)

      character(300), intent(in) :: file

      restart_file = file

   end subroutine init_save

   !****************************************************************************
   !                                      Title
   !****************************************************************************
   subroutine init_dt(delt0, istatus)

#ifdef NETCDF
      use mp, only: proc0, broadcast
      use file_utils, only: error_unit
#endif
      implicit none
      real, intent(in out) :: delt0
      integer, intent(out) :: istatus
#ifdef NETCDF
      character(306) :: path_netcdf_file_per_proc

      if (proc0) then

         if (.not. initialised_restart_module) then

#ifdef NETCDF_PARALLEL
            if (save_many) then
#endif
               path_netcdf_file_per_proc = trim(trim(restart_file)//'.0')
#ifdef NETCDF_PARALLEL
            else
               path_netcdf_file_per_proc = trim(trim(restart_file))
            end if
#endif

            istatus = nf90_open(path_netcdf_file_per_proc, NF90_NOWRITE, ncid)
            if (istatus /= NF90_NOERR) call netcdf_error(istatus, file=path_netcdf_file_per_proc)

            istatus = nf90_inq_varid(ncid, "delt0", delt0id)
            if (istatus /= NF90_NOERR) call netcdf_error(istatus, var='delt0')
         end if

         istatus = nf90_get_var(ncid, delt0id, delt0)

         if (istatus /= NF90_NOERR) then
            call netcdf_error(istatus, ncid, delt0id, message=' in init_dt')
            delt0 = -1.
         end if

         if (.not. initialised_restart_module) istatus = nf90_close(ncid)
      end if

      call broadcast(istatus)
      call broadcast(delt0)

#endif

   end subroutine init_dt

   !****************************************************************************
   !                                      Title
   !****************************************************************************
   subroutine init_tstart(tstart, istep0, istatus)

#ifdef NETCDF
      use mp, only: proc0, broadcast, mp_abort
      use file_utils, only: error_unit
#endif

      use parameters_numerical, only: nstep
      
      implicit none
      real, intent(in out) :: tstart
      integer, intent(out) :: istep0
      integer, intent(out) :: istatus
#ifdef NETCDF
      character(306) :: path_netcdf_file_per_proc

      if (proc0) then
#ifdef NETCDF_PARALLEL
         if (save_many) then
#endif
            path_netcdf_file_per_proc = trim(trim(restart_file)//'.0')
#ifdef NETCDF_PARALLEL
         else
            path_netcdf_file_per_proc = trim(trim(restart_file))
         end if
#endif

         if (.not. initialised_restart_module) then

            istatus = nf90_open(path_netcdf_file_per_proc, NF90_NOWRITE, ncid)
            if (istatus /= NF90_NOERR) call netcdf_error(istatus, file=path_netcdf_file_per_proc)
         end if

         istatus = nf90_inq_varid(ncid, "t0", t0id)
         if (istatus /= NF90_NOERR) call netcdf_error(istatus, var='t0')

         istatus = nf90_get_var(ncid, t0id, tstart)
         if (istatus /= NF90_NOERR) then
            call netcdf_error(istatus, ncid, t0id, message=' in init_tstart')
            tstart = -1.
         end if

         istatus = nf90_inq_varid(ncid, "istep0", istep0id)
         if (istatus /= NF90_NOERR) call netcdf_error(istatus, var='istep0')

         istatus = nf90_get_var(ncid, istep0id, istep0)
         if (istatus /= NF90_NOERR) then
            call netcdf_error(istatus, ncid, istep0id, message=' in init_tstart')
            istep0 = -1
         end if

         if (.not. initialised_restart_module) istatus = nf90_close(ncid)

      end if
      
      ! Sanity checks
      if (nstep > 0) then
         if (nstep <= istep0) then
            call mp_abort('Restarted simulation has nstep < istep. Aborting.')
         end if
      end if 
      
      call broadcast(istatus)
      call broadcast(istep0)
      call broadcast(tstart)

#endif

   end subroutine init_tstart
   
   
   !****************************************************************************
   !                        Abort if an error was encountered                   
   !****************************************************************************
   subroutine process_nf90_error(error_message, istatus)
   
      use mp, only: mp_abort, proc0
      use file_units, only: unit_error_file
     
      integer, intent(in) :: istatus
      character(*), intent(in) :: error_message
      
      !-------------------------------------------------------------------------
      
      if (.not. proc0) return
      write (unit_error_file, *) error_message, nf90_strerror(istatus)
      write (*, *) ' '; write (*, *) error_message, nf90_strerror(istatus)
      call mp_abort('Error while writing the netcdf file to restart a simulation. Aborting.')
      
   end subroutine process_nf90_error

end module save_stella_for_restart
