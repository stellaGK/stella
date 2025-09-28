# include "define.inc"
!###############################################################################
!                   SAVE DISTRIBUTION FUNCTION FOR A RESTART                    
!###############################################################################
! This module saves the distribution function g(mu,vpa,ikxkyzs) to a netcdf
! file, as well as the dimensions (kx,ky,z,nu,vpa,s,tube). This allows us to
! restart the simulation in the future, to continue the time evolution.
! 
! Note that we also need to save the <shift_state> for Flow shear, otherwise
! errors will pop up. For radial variation simulations or for simulations 
! including sources, we also need to save <int_krook>, <int_proj>, <g_krook>,
! <g_proj> and <phi_proj>, in order to restart a simulation.
! 
! WARNING: This module has not been cleaned completely, since parallel netcdf 
! does not seem to work. Therefore, one should always run with <save_many> = .true.
! until this module has been debugged and fixed.
!###############################################################################
module save_stella_for_restart_parallel_netcdf

   ! Import mpi
   use mp, only: mp_comm
   use mp, only: mp_info

   ! Import parallel netcdf modules
   ! If using netcdf version 4.1.2 or older replace NF90_MPIIO with NF90_CLOBBER
#ifdef NETCDF_PARALLEL
   use netcdf, only: NF90_HDF5, NF90_MPIIO
   use netcdf, only: nf90_var_par_access, NF90_COLLECTIVE
   use netcdf, only: nf90_put_att, NF90_GLOBAL, nf90_get_att
   use netcdf, only: NF90_NOWRITE, NF90_NOERR
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
#ifdef NETCDF_PARALLEL
   public :: save_stella_data_for_restart_to_a_single_file
   public :: read_stella_data_for_restart_from_single_file
#endif

   private

   ! Initialise local arrays
# ifdef NETCDF_PARALLEL
   real, allocatable, dimension(:, :, :) :: tmpr, tmpi
   real, allocatable, dimension(:, :, :, :) :: ktmpr, ktmpi
   real, allocatable, dimension(:, :, :, :)   :: ptmpr, ptmpi
   real, allocatable, dimension(:, :, :)   :: pptmpr, pptmpi
   integer(kind_nf) :: ncid, zedid, vpaid, gloid, gvmuloid, kyid, kxid, muid, tubeid
   integer(kind_nf) :: krookr_id, krooki_id, projr_id, proji_id
   integer(kind_nf) :: phiprojr_id, phiproji_id
   integer(kind_nf) :: t0id, gr_id, gi_id, delt0id, istep0id
   integer(kind_nf) :: intkrook_id, intproj_id;
   integer(kind_nf) :: shift_id
   logical :: initialised_restart_module = .false.
#endif

contains


!###############################################################################
!############################### SAVE DISTRIBUTION #############################
!###############################################################################

#ifdef NETCDF_PARALLEL
   !****************************************************************************
   !                 Save distribution function to a netcdf file                
   !****************************************************************************
   subroutine save_stella_data_for_restart_to_a_single_file(istep0, t0, delt0, istatus, exit_in, restart_file)
      
      ! Must include kxkyz_layout_type here to avoid obscure bomb while compiling
      ! diagnostics.f90 (which uses this module) with the Compaq F90 compiler:
      use parallelisation_layouts, only: kxkyz_lo, vmu_lo
      use common_types, only: kxkyz_layout_type
      use mp, only: iproc, barrier
      use mp, only: proc0
      
      ! Grids
      use grids_z, only: nztot, nzgrid, ntubes
      use grids_velocity, only: nvpa, nmu
      use grids_kxky, only: naky, nakx
      
      ! Radial variation and sources
      use gk_sources, only: include_qn_source
      use gk_sources, only: source_option_krook, source_option_projection
      use gk_sources, only: source_option_switch, int_krook, int_proj
      use arrays_distribution_function, only: g_krook, g_proj
      use arrays_fields, only: phi_proj
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
      character(300), intent(in out) :: restart_file
      
      ! Local variables
      character(306) :: path_netcdf_file
      integer :: i, n_elements, nvmulo_elements, tmpunit
      integer :: total_elements, total_vmulo_elements
      logical :: has_vmulo
      integer, dimension(3) :: start_pos, counts
      integer, dimension(4) :: start_pos_krook, counts_krook
      logical :: exit

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
      total_elements = kxkyz_lo%ulim_world + 1
      nvmulo_elements = vmu_lo%ulim_proc - vmu_lo%llim_proc + 1
      total_vmulo_elements = vmu_lo%ulim_world + 1
      if (n_elements <= 0) return

      ! Create netcdf file for each vmulo element
      has_vmulo = .true.

      !-------------------------------------------------------------------------
      !                      Initialise the restart module                      
      !-------------------------------------------------------------------------
      if (.not. initialised_restart_module) then

         ! Only initialise once
         initialised_restart_module = .true.
         
         ! Path of the netcdf file
         path_netcdf_file = trim(restart_file)

         ! Create a single netcdf file using the parallel netcdf library
         call barrier
         tmpunit = 983
         if (iproc == 0) open (unit=tmpunit, file=path_netcdf_file)
         if (iproc == 0) close (unit=tmpunit, status='delete')
         call barrier
         istatus = nf90_create(path_netcdf_file, IOR(NF90_HDF5, NF90_MPIIO), ncid, comm=mp_comm, info=mp_info)
         if (istatus /= NF90_NOERR) call process_nf90_error("nf90_create error: ", istatus)

         ! Save the dimensions
         if (n_elements > 0) then
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
            istatus = nf90_def_dim(ncid, "glo", total_elements, gloid)
            if (istatus /= NF90_NOERR) call process_nf90_error("nf90_def_dim glo error: ", istatus)
            istatus = nf90_def_dim(ncid, "gvmulo", total_vmulo_elements, gvmuloid)
            if (istatus /= NF90_NOERR) call process_nf90_error("nf90_def_dim gvmulo error: ", istatus)
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
            if (istatus /= NF90_NOERR) call process_nf90_error("nf90_def_var g real error: ", istatus)
            istatus = nf90_def_var(ncid, "gi", netcdf_real, (/vpaid, muid, gloid/), gi_id)
            if (istatus /= NF90_NOERR) call process_nf90_error("nf90_def_var g imag error: ", istatus)

            ! Flow shear
            istatus = nf90_def_var(ncid, "shiftstate", netcdf_real, (/kyid/), shift_id)
            if (istatus /= NF90_NOERR) call process_nf90_error("nf90_def_var shiftstate error: ", istatus)
 
            ! Radial variation and Sources
            if (source_option_switch == source_option_krook .and. has_vmulo) then
               istatus = nf90_def_var(ncid, "intkrook", netcdf_real, intkrook_id)
               if (istatus /= NF90_NOERR) call process_nf90_error("nf90_def_var intkrook error: ", istatus)
               istatus = nf90_def_var(ncid, "krookr", netcdf_real, (/kxid, zedid, tubeid, gvmuloid/), krookr_id)
               if (istatus /= NF90_NOERR) call process_nf90_error("nf90_def_var krookr error: ", istatus)
               istatus = nf90_def_var(ncid, "krooki", netcdf_real, (/kxid, zedid, tubeid, gvmuloid/), krooki_id)
               if (istatus /= NF90_NOERR) call process_nf90_error("nf90_def_var krooki error: ", istatus)
            end if

            ! Radial variation and Sources
            if (include_qn_source .and. iproc == 0) then
               istatus = nf90_def_var(ncid, "phiprojr", netcdf_real, (/kxid, zedid, tubeid/), phiprojr_id)
               if (istatus /= NF90_NOERR) call process_nf90_error("nf90_def_var phiprojr error: ", istatus)
               istatus = nf90_def_var(ncid, "phiproji", netcdf_real, (/kxid, zedid, tubeid/), phiproji_id)
               if (istatus /= NF90_NOERR) call process_nf90_error("nf90_def_var phiproji error: ", istatus)
            end if

            ! Radial variation and Sources
            if (source_option_switch == source_option_projection .and. has_vmulo) then
               istatus = nf90_def_var(ncid, "intproj", netcdf_real, intproj_id)
               if (istatus /= NF90_NOERR) call process_nf90_error("nf90_def_var intproj error: ", istatus)
               istatus = nf90_def_var(ncid, "projr", netcdf_real, (/kxid, zedid, tubeid, gvmuloid/), projr_id)
               if (istatus /= NF90_NOERR) call process_nf90_error("nf90_def_var projr error: ", istatus)
               istatus = nf90_def_var(ncid, "proji", netcdf_real, (/kxid, zedid, tubeid, gvmuloid/), proji_id)
               if (istatus /= NF90_NOERR) call process_nf90_error("nf90_def_var proji error: ", istatus)
            end if
         end if

         ! Finished defining all dimensions and variables in the netcdf file
         istatus = nf90_enddef(ncid)
         if (istatus /= NF90_NOERR) call process_nf90_error("nf90_enddef error: ", istatus)
         
      end if
      
      !-------------------------------------------------------------------------
      !                          Write the netcdf file                          
      !-------------------------------------------------------------------------

      ! Only write on the first processor
      if (iproc == 0) then

         ! Write the actual time, istep and delt variables to the netcdf file
         istatus = nf90_put_var(ncid, delt0id, delt0)
         if (istatus /= NF90_NOERR) call process_nf90_error("nf90_put_var delt0 error: ", istatus)
         istatus = nf90_put_var(ncid, t0id, t0)
         if (istatus /= NF90_NOERR) call process_nf90_error("nf90_put_var t0 error: ", istatus)
         istatus = nf90_put_var(ncid, istep0id, istep0)
         if (istatus /= NF90_NOERR) call process_nf90_error("nf90_put_var istep0 error: ", istatus)
      end if

      ! Save the distribution function to the netcdf file
      if (n_elements > 0) then
      
         ! Check whether the variables have been initialised
         istatus = nf90_var_par_access(ncid, gr_id, NF90_COLLECTIVE)
         if (istatus /= NF90_NOERR) call netcdf_error(istatus, ncid, gr_id)
         istatus = nf90_var_par_access(ncid, gi_id, NF90_COLLECTIVE)
         if (istatus /= NF90_NOERR) call netcdf_error(istatus, ncid, gr_id)
         
         ! Dimensions on the parallelised grid
         start_pos = (/1, 1, kxkyz_lo%llim_proc + 1/)
         counts = (/nvpa, nmu, n_elements/)

         ! Allocate temporary array
         if (.not. allocated(tmpr)) allocate (tmpr(nvpa, nmu, kxkyz_lo%llim_proc:kxkyz_lo%ulim_alloc))

         ! Save the real part of the distribution function
         tmpr = real(gvmu)
         istatus = nf90_put_var(ncid, gr_id, tmpr, start=start_pos, count=counts)
         if (istatus /= NF90_NOERR) call netcdf_error(istatus, ncid, gr_id)

         ! Save the imaginary part of the distribution function
         tmpr = aimag(gvmu)
         istatus = nf90_put_var(ncid, gi_id, tmpr, start=start_pos, count=counts)
         if (istatus /= NF90_NOERR) call netcdf_error(istatus, ncid, gi_id)

         ! Flow shear
         istatus = nf90_put_var(ncid, shift_id, shift_state)
         if (istatus /= NF90_NOERR) call netcdf_error(istatus, ncid, shift_id)

         ! Radial variation and sources
         if (source_option_switch == source_option_krook .and. has_vmulo) then
            if (.not. allocated(ktmpr)) allocate (ktmpr(nakx, -nzgrid:nzgrid, ntubes, vmu_lo%llim_proc:vmu_lo%ulim_alloc))
            if (.not. allocated(ktmpi)) allocate (ktmpi(nakx, -nzgrid:nzgrid, ntubes, vmu_lo%llim_proc:vmu_lo%ulim_alloc))
            if (iproc == 0) then
               istatus = nf90_put_var(ncid, intkrook_id, int_krook)
               if (istatus /= NF90_NOERR) call process_nf90_error("nf90_put_var int_krook error: ", istatus)
            end if
            ktmpr = real(g_krook)
            ktmpi = aimag(g_krook)
            istatus = nf90_var_par_access(ncid, krookr_id, NF90_COLLECTIVE)
            istatus = nf90_var_par_access(ncid, krooki_id, NF90_COLLECTIVE)
            start_pos_krook = (/1, 1, 1, vmu_lo%llim_proc + 1/)
            counts_krook = (/nakx, nztot, ntubes, nvmulo_elements/)
            istatus = nf90_put_var(ncid, krookr_id, ktmpr, start=start_pos_krook, count=counts_krook)
            if (istatus /= NF90_NOERR) call netcdf_error(istatus, ncid, krookr_id)
            istatus = nf90_put_var(ncid, krooki_id, ktmpi, start=start_pos_krook, count=counts_krook)
            if (istatus /= NF90_NOERR) call netcdf_error(istatus, ncid, krooki_id)
         end if

         ! Radial variation and sources
         if (source_option_switch == source_option_projection .and. has_vmulo) then
            if (.not. allocated(ptmpr)) allocate (ptmpr(nakx, -nzgrid:nzgrid, ntubes, vmu_lo%llim_proc:vmu_lo%ulim_alloc))
            if (.not. allocated(ptmpi)) allocate (ptmpi(nakx, -nzgrid:nzgrid, ntubes, vmu_lo%llim_proc:vmu_lo%ulim_alloc))
            if (iproc == 0) then
               istatus = nf90_put_var(ncid, intproj_id, int_proj)
               if (istatus /= NF90_NOERR) call process_nf90_error("nf90_put_var int_proj error: ", istatus)
            end if
            ptmpr = real(g_proj)
            ptmpi = aimag(g_proj)
            istatus = nf90_var_par_access(ncid, projr_id, NF90_COLLECTIVE)
            istatus = nf90_var_par_access(ncid, proji_id, NF90_COLLECTIVE)
            start_pos_krook = (/1, 1, 1, vmu_lo%llim_proc + 1/)
            counts_krook = (/nakx, nztot, ntubes, nvmulo_elements/)
            istatus = nf90_put_var(ncid, projr_id, ptmpr, start=start_pos_krook, count=counts_krook)
            if (istatus /= NF90_NOERR) call netcdf_error(istatus, ncid, projr_id)
            istatus = nf90_put_var(ncid, proji_id, ptmpi, start=start_pos_krook, count=counts_krook)
            if (istatus /= NF90_NOERR) call netcdf_error(istatus, ncid, proji_id)
         end if

         ! Radial variation and sources
         if (include_qn_source .and. iproc == 0) then
            if (.not. allocated(pptmpr)) allocate (pptmpr(nakx, -nzgrid:nzgrid, ntubes))
            if (.not. allocated(pptmpi)) allocate (pptmpi(nakx, -nzgrid:nzgrid, ntubes))
            pptmpr = real(phi_proj)
            pptmpi = aimag(phi_proj)
            istatus = nf90_var_par_access(ncid, phiprojr_id, NF90_COLLECTIVE)
            istatus = nf90_var_par_access(ncid, phiproji_id, NF90_COLLECTIVE)
            start_pos = (/1, 1, 1/)
            counts = (/nakx, nztot, ntubes/)
            istatus = nf90_put_var(ncid, phiprojr_id, ktmpr, start=start_pos, count=counts)
            if (istatus /= NF90_NOERR) call netcdf_error(istatus, ncid, phiprojr_id)
            istatus = nf90_put_var(ncid, phiproji_id, pptmpi, start=start_pos, count=counts)
            if (istatus /= NF90_NOERR) call netcdf_error(istatus, ncid, phiproji_id)
         end if
         
      end if

      ! If we make a clean exit of stella, then close the netcdf file
      if (exit) then
         i = nf90_close(ncid)
         if (i /= NF90_NOERR) call netcdf_error(istatus, message='nf90_close error')
         
      ! Otherwise, sync the netcdf file
      else
         i = nf90_sync(ncid)
         if (i /= NF90_NOERR) call netcdf_error(istatus, message='nf90_sync error')
      end if

      ! Deallocate temporary arrays
      if (allocated(tmpr)) deallocate (tmpr)
      if (allocated(tmpi)) deallocate (tmpi)
      if (allocated(ptmpr)) deallocate (ptmpr)
      if (allocated(ptmpi)) deallocate (ptmpi)
      if (allocated(ktmpr)) deallocate (ktmpr)
      if (allocated(ktmpi)) deallocate (ktmpi)
      if (allocated(pptmpr)) deallocate (pptmpr)
      if (allocated(pptmpi)) deallocate (pptmpi)

   end subroutine save_stella_data_for_restart_to_a_single_file
#endif

!###############################################################################
!############################### READ DISTRIBUTION #############################
!###############################################################################

#ifdef NETCDF_PARALLEL
   !****************************************************************************
   !            Read distribution function from a single netcdf file            
   !****************************************************************************
   subroutine read_stella_data_for_restart_from_single_file(g, scale, istatus, restart_file)
      
      ! Parallelisation
      use mp, only: iproc, broadcast, proc0
      use parallelisation_layouts, only: kxkyz_lo, vmu_lo
      
      ! Grids
      use grids_z, only: nztot, nzgrid, ntubes
      use grids_velocity, only: nvpa, nmu
      use grids_kxky, only: naky, nakx
      
      ! Radial variation and sources
      use gk_sources, only: source_option_krook, source_option_projection
      use gk_sources, only: source_option_switch, int_krook, int_proj
      use gk_sources, only: include_qn_source
      use arrays_distribution_function, only: g_krook, g_proj
      use arrays_fields, only: phi_proj
      
      ! Error files
      use file_units, only: unit_error_file
      
      ! Flow shear
      use arrays, only: shift_state

      implicit none

      ! Arguments
      complex, dimension(:, :, kxkyz_lo%llim_proc:), intent(out) :: g
      real, intent(in) :: scale
      integer, intent(out) :: istatus
      character(300), intent(in out) :: restart_file
   
      ! Variables
      integer, dimension(3) :: counts, start_pos
      integer, dimension(4) :: counts_krook, start_pos_krook
      character(306) :: path_netcdf_file_per_proc
      integer :: i, n_elements, nvmulo_elements
      logical :: has_vmulo
      
      real, allocatable, dimension(:, :, :) :: tmpr, tmpi
      real, allocatable, dimension(:, :, :, :) :: ktmpr, ktmpi
      real, allocatable, dimension(:, :, :, :) :: ptmpr, ptmpi
      real, allocatable, dimension(:, :, :) :: pptmpr, pptmpi

      n_elements = kxkyz_lo%ulim_proc - kxkyz_lo%llim_proc + 1
      nvmulo_elements = vmu_lo%ulim_proc - vmu_lo%llim_proc + 1

      if (n_elements <= 0) return

      has_vmulo = .true.

      if (.not. initialised_restart_module) then
!       initialised_restart_module = .true.
         path_netcdf_file_per_proc = trim(restart_file)

! If using netcdf version 4.1.2 deleted NF90_MPIIO and the associated IOR
         istatus = nf90_open(path_netcdf_file_per_proc, IOR(NF90_NOWRITE, NF90_MPIIO), ncid, comm=mp_comm, info=mp_info)
         if (istatus /= NF90_NOERR) then
            call netcdf_error(istatus, file=path_netcdf_file_per_proc, abort=.true.)
         end if

         ! check precision
         if (netcdf_real == 0) netcdf_real = get_netcdf_code_precision()
         call check_netcdf_file_precision(ncid)

         istatus = nf90_inq_dimid(ncid, "tube", tubeid)
         if (istatus /= NF90_NOERR) call netcdf_error(istatus, dim='tube')

         istatus = nf90_inq_dimid(ncid, "zed", zedid)
         if (istatus /= NF90_NOERR) call netcdf_error(istatus, dim='zed')

         istatus = nf90_inq_dimid(ncid, "aky", kyid)
         if (istatus /= NF90_NOERR) call netcdf_error(istatus, dim='aky')

         istatus = nf90_inq_dimid(ncid, "akx", kxid)
         if (istatus /= NF90_NOERR) call netcdf_error(istatus, dim='akx')

         istatus = nf90_inq_dimid(ncid, "glo", gloid)
         if (istatus /= NF90_NOERR) call netcdf_error(istatus, dim='glo')

         if (has_vmulo) then
            istatus = nf90_inq_dimid(ncid, "gvmulo", gvmuloid)
            if (istatus /= NF90_NOERR) call netcdf_error(istatus, dim='gvmulo')
         end if

         istatus = nf90_inquire_dimension(ncid, tubeid, len=i)
         if (istatus /= NF90_NOERR) call netcdf_error(istatus, ncid, dimid=tubeid)
         if (i /= ntubes) write (*, *) 'Restart error: ntubes=? ', i, ' : ', ntubes, ' : ', iproc

         istatus = nf90_inquire_dimension(ncid, zedid, len=i)
         if (istatus /= NF90_NOERR) call netcdf_error(istatus, ncid, dimid=zedid)
         if (i /= 2 * nzgrid + 1) write (*, *) 'Restart error: nzgrid=? ', i, ' : ', nzgrid, ' : ', iproc

         istatus = nf90_inquire_dimension(ncid, kyid, len=i)
         if (istatus /= NF90_NOERR) call netcdf_error(istatus, ncid, dimid=kyid)
         if (i /= naky) write (*, *) 'Restart error: naky=? ', i, ' : ', naky, ' : ', iproc

         istatus = nf90_inquire_dimension(ncid, kxid, len=i)
         if (istatus /= NF90_NOERR) call netcdf_error(istatus, ncid, dimid=kxid)
         if (i /= nakx) write (*, *) 'Restart error: nakx=? ', i, ' : ', nakx, ' : ', iproc

         istatus = nf90_inquire_dimension(ncid, gloid, len=i)
         if (istatus /= NF90_NOERR) call netcdf_error(istatus, ncid, dimid=gloid)

         if (i /= kxkyz_lo%ulim_world + 1) write (*, *) 'Restart error: glo=? ', i, ' : ', iproc

         if (has_vmulo) then
            istatus = nf90_inquire_dimension(ncid, gvmuloid, len=i)
            if (istatus /= NF90_NOERR) call netcdf_error(istatus, ncid, dimid=gvmuloid)
            if (i /= vmu_lo%ulim_world + 1) write (*, *) 'Restart error: gvmulo=? ', i, ' : ', iproc
         end if

         if (source_option_switch == source_option_krook .and. has_vmulo) then
            istatus = nf90_inq_varid(ncid, "intkrook", intkrook_id)
            if (istatus /= NF90_NOERR) call netcdf_error(istatus, var='intkrook')

            istatus = nf90_inq_varid(ncid, "krookr", krookr_id)
            if (istatus /= NF90_NOERR) call netcdf_error(istatus, var='krookr')

            istatus = nf90_inq_varid(ncid, "krooki", krooki_id)
            if (istatus /= NF90_NOERR) call netcdf_error(istatus, var='krooki')

         end if

         if (source_option_switch == source_option_projection .and. has_vmulo) then
            istatus = nf90_inq_varid(ncid, "intproj", intproj_id)
            if (istatus /= NF90_NOERR) call netcdf_error(istatus, var='intproj')

            istatus = nf90_inq_varid(ncid, "projr", projr_id)
            if (istatus /= NF90_NOERR) call netcdf_error(istatus, var='projr')

            istatus = nf90_inq_varid(ncid, "proji", proji_id)
            if (istatus /= NF90_NOERR) call netcdf_error(istatus, var='proji')
         end if

         if (include_qn_source .and. iproc == 0) then
            istatus = nf90_inq_varid(ncid, "phiprojr", phiprojr_id)
            if (istatus /= NF90_NOERR) call netcdf_error(istatus, var='phiprojr')

            istatus = nf90_inq_varid(ncid, "phiproji", phiproji_id)
            if (istatus /= NF90_NOERR) call netcdf_error(istatus, var='phiproji')
         end if

         istatus = nf90_inq_varid(ncid, "shiftstate", shift_id)
         if (istatus /= NF90_NOERR) call netcdf_error(istatus, var='shiftstate')

         istatus = nf90_inq_varid(ncid, "gr", gr_id)
         if (istatus /= NF90_NOERR) call netcdf_error(istatus, var='gr')

         istatus = nf90_inq_varid(ncid, "gi", gi_id)
         if (istatus /= NF90_NOERR) call netcdf_error(istatus, var='gi')
      end if

      if (.not. allocated(tmpr)) allocate (tmpr(nvpa, nmu, kxkyz_lo%llim_proc:kxkyz_lo%ulim_alloc))
      if (.not. allocated(tmpi)) allocate (tmpi(nvpa, nmu, kxkyz_lo%llim_proc:kxkyz_lo%ulim_alloc))

      tmpr = 0.; tmpi = 0.
      start_pos = (/1, 1, kxkyz_lo%llim_proc + 1/)
      counts = (/nvpa, nmu, n_elements/)
      istatus = nf90_get_var(ncid, gr_id, tmpr, start=start_pos, count=counts)

      if (istatus /= NF90_NOERR) call netcdf_error(istatus, ncid, gr_id)

      istatus = nf90_get_var(ncid, gi_id, tmpi, start=start_pos_krook, count=counts_krook)

      if (istatus /= NF90_NOERR) call netcdf_error(istatus, ncid, gi_id)

      g = cmplx(tmpr, tmpi)

      if (source_option_switch == source_option_krook .and. has_vmulo) then
         if (.not. allocated(ktmpr)) &
            allocate (ktmpr(nakx, -nzgrid:nzgrid, ntubes, vmu_lo%llim_proc:vmu_lo%ulim_alloc))
         if (.not. allocated(ktmpi)) &
            allocate (ktmpi(nakx, -nzgrid:nzgrid, ntubes, vmu_lo%llim_proc:vmu_lo%ulim_alloc))

         istatus = nf90_get_var(ncid, intkrook_id, int_krook)
         if (istatus /= NF90_NOERR) call netcdf_error(istatus, ncid, intkrook_id)

         ktmpr = 0.; ktmpi = 0.
         start_pos_krook = (/1, 1, 1, vmu_lo%llim_proc + 1/)
         counts_krook = (/nakx, nztot, ntubes, nvmulo_elements/)
         istatus = nf90_get_var(ncid, krookr_id, ktmpr, start=start_pos_krook, count=counts_krook)
         if (istatus /= NF90_NOERR) call netcdf_error(istatus, ncid, krookr_id)
         istatus = nf90_get_var(ncid, krooki_id, ktmpi, start=start_pos_krook, count=counts_krook)
         if (istatus /= NF90_NOERR) call netcdf_error(istatus, ncid, krooki_id)

         g_krook = cmplx(ktmpr, ktmpi)

      end if

      if (source_option_switch == source_option_projection .and. has_vmulo) then
         if (.not. allocated(ptmpr)) &
            allocate (ptmpr(nakx, -nzgrid:nzgrid, ntubes, vmu_lo%llim_proc:vmu_lo%ulim_alloc))
         if (.not. allocated(ptmpi)) &
            allocate (ptmpi(nakx, -nzgrid:nzgrid, ntubes, vmu_lo%llim_proc:vmu_lo%ulim_alloc))

         istatus = nf90_get_var(ncid, intproj_id, int_proj)
         if (istatus /= NF90_NOERR) call netcdf_error(istatus, ncid, intproj_id)

         ptmpr = 0.; ptmpi = 0.
         start_pos_krook = (/1, 1, 1, vmu_lo%llim_proc + 1/)
         counts_krook = (/nakx, nztot, ntubes, nvmulo_elements/)
         istatus = nf90_get_var(ncid, projr_id, ptmpr, start=start_pos_krook, count=counts_krook)

         if (istatus /= NF90_NOERR) call netcdf_error(istatus, ncid, projr_id)

         istatus = nf90_get_var(ncid, proji_id, ptmpi, start=start_pos_krook, count=counts_krook)

         if (istatus /= NF90_NOERR) call netcdf_error(istatus, ncid, proji_id)

         g_proj = cmplx(ptmpr, ptmpi)

      end if

      if (include_qn_source .and. iproc == 0) then
         if (.not. allocated(pptmpr)) allocate (pptmpr(nakx, -nzgrid:nzgrid, ntubes))
         if (.not. allocated(pptmpi)) allocate (pptmpi(nakx, -nzgrid:nzgrid, ntubes))

         pptmpr = 0.; pptmpi = 0.
         start_pos = (/1, 1, 1/)
         counts = (/nakx, nztot, ntubes/)
         istatus = nf90_get_var(ncid, phiprojr_id, pptmpr, start=start_pos, count=counts)

         if (istatus /= NF90_NOERR) call netcdf_error(istatus, ncid, phiprojr_id)

         istatus = nf90_get_var(ncid, phiproji_id, pptmpi, start=start_pos, count=counts)

         if (istatus /= NF90_NOERR) call netcdf_error(istatus, ncid, phiproji_id)

         phi_proj = cmplx(pptmpr, pptmpi)

      end if

      if (.not. allocated(shift_state)) allocate (shift_state(naky))
      istatus = nf90_get_var(ncid, shift_id, shift_state)
      if (istatus /= NF90_NOERR) call netcdf_error(istatus, ncid, shift_id)

      if (scale > 0.) then
         g = g * scale
         if (source_option_switch == source_option_krook) g_krook = g_krook * scale
         if (source_option_switch == source_option_projection) g_proj = g_proj * scale
      end if

      ! RN 2008/05/23: this was commented out. why? HJL 2013/05/15 Because it stops future writing to the file
!    istatus = nf90_close (ncid)
      if (istatus /= NF90_NOERR) then
         write (unit_error_file, *) "nf90_close error: ", nf90_strerror(istatus), ' ', iproc
      end if

      if (allocated(tmpr)) deallocate (tmpr)
      if (allocated(tmpi)) deallocate (tmpi)
      if (allocated(ptmpr)) deallocate (ptmpr)
      if (allocated(ptmpi)) deallocate (ptmpi)
      if (allocated(ktmpr)) deallocate (ktmpr)
      if (allocated(ktmpi)) deallocate (ktmpi)
      if (allocated(pptmpr)) deallocate (pptmpr)
      if (allocated(pptmpi)) deallocate (pptmpi)

      if (include_qn_source) call broadcast(phi_proj)

   end subroutine read_stella_data_for_restart_from_single_file
#endif

#ifdef NETCDF_PARALLEL
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
#endif

end module save_stella_for_restart_parallel_netcdf
