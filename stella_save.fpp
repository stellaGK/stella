# include "define.inc"

module stella_save

  use mp, only: mp_comm, mp_info

# ifdef NETCDF
!  use netcdf, only: NF90_FLOAT, NF90_DOUBLE
# ifdef NETCDF_PARALLEL
! If using netcdf version 4.1.2 or older delete NF90_MPIIO
  use netcdf, only: NF90_HDF5,NF90_MPIIO
  use netcdf, only: nf90_var_par_access, NF90_COLLECTIVE
  use netcdf, only: nf90_put_att, NF90_GLOBAL, nf90_get_att
# endif
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
# endif

  implicit none

  public :: stella_restore, stella_save_for_restart
  public :: read_many, save_many
  public :: init_save, init_dt, init_tstart, finish_save

!# ifdef NETCDF
!  public :: netcdf_real, kind_nf, get_netcdf_code_precision, netcdf_error
!# endif

  interface stella_restore
     module procedure stella_restore_many
  end interface

  logical :: read_many, save_many ! Read and write single or multiple restart files

  private
  character (300), save :: restart_file

# ifdef NETCDF
  real, allocatable, dimension (:,:,:) :: tmpr, tmpi
  real, allocatable, dimension (:,:,:) :: ktmpr, ktmpi
  real, allocatable, dimension (:,:,:)   :: ptmpr, ptmpi
  real, allocatable, dimension (:,:,:,:) :: ftmpr, ftmpi
  integer (kind_nf) :: ncid, zedid, vpaid, gloid, gvmuloid, kyid, kxid, muid, tubeid
  integer (kind_nf) :: phir_id, phii_id, aparr_id, apari_id, bparr_id, bpari_id ! Bob: added here
  integer (kind_nf) :: krookr_id, krooki_id, projr_id, proji_id
  integer (kind_nf) :: t0id, gr_id, gi_id, delt0id, istep0id
  integer (kind_nf) :: intkrook_id, intproj_id;
  integer (kind_nf) :: shift_id

  logical :: initialized = .false.
# endif

contains

!!----------------------------------------------------------------------!!
!!----------------------------------------------------------------------!!
!!--Save----------------------------------------------------------------!!
!!----------------------------------------------------------------------!!
!!----------------------------------------------------------------------!!

  subroutine stella_save_for_restart &
       (g, istep0, t0, delt0, istatus, fphi, fapar, fbpar, exit_in, fileopt)

# ifdef NETCDF
    use fields_arrays, only: phi, apar, bpar, shift_state
    use dist_fn_arrays, only: g_krook, g_proj
    use kt_grids, only: naky, nakx
# else
    use mp, only: proc0
# endif
    use mp, only: iproc, barrier
    use zgrid, only: nzgrid, ntubes
    ! Must include kxkyz_layout_type here to avoid obscure bomb while compiling
    ! stella_diagnostics.f90 (which uses this module) with the Compaq F90 compiler:
    use stella_layouts, only: kxkyz_lo, xyzs_layout, vms_layout, vmu_lo
    use common_types, only: kxkyz_layout_type
    use file_utils, only: error_unit
    use vpamu_grids, only: nvpa, nmu
    use dissipation, only: include_krook_operator, int_krook
    use dissipation, only: remove_zero_projection, int_proj
    use physics_flags, only: prp_shear_enabled

    implicit none

    complex, dimension (:,:,kxkyz_lo%llim_proc:), intent (in) :: g
    real, intent (in) :: t0, delt0
    real, intent (in) :: fphi, fapar, fbpar
    integer, intent (in) :: istep0
    integer, intent (out) :: istatus
    logical, intent (in), optional :: exit_in
    character (20), INTENT (in), optional :: fileopt
# ifdef NETCDF
    character (306) :: file_proc
    character (10) :: suffix
    integer :: i, n_elements, nvmulo_elements, ierr
    integer :: total_elements, total_vmulo_elements
# ifdef NETCDF_PARALLEL
    integer, dimension(3) :: start_pos, counts
# endif
    logical :: exit

!*********-----------------------_**********************

    istatus = 0
    if (present(exit_in)) then
       exit = exit_in
    else
       exit = .false.
    end if

!    if (proc0) then
!      write (*,*) "Starting save_for_restart in ", restart_file
!      write (*,*) "List restart files"
!      call system("echo 'start' >> filelist.txt; ls nc/* >> filelist.txt;  ")
!    end if

    n_elements = kxkyz_lo%ulim_proc-kxkyz_lo%llim_proc+1
    total_elements = kxkyz_lo%ulim_world+1

    nvmulo_elements = vmu_lo%ulim_proc-vmu_lo%llim_proc+1
    total_vmulo_elements = vmu_lo%ulim_world+1

    if (n_elements <= 0) return

    if (.not.initialized) then

       initialized = .true.

       file_proc = trim(restart_file)

!CMR, 5/4/2011: Add optional piece of filename
       IF (PRESENT(fileopt)) THEN
          file_proc=trim(file_proc)//trim(fileopt)
       END IF
!CMRend

!</HL>  The NETCDF_PARALLEL directives include code for parallel
!       netcdf using HDF5 to write the output to a single restart file
!       The read_many flag allows the old style multiple file output
# ifdef NETCDF_PARALLEL
       if(save_many) then
# endif
          WRITE (suffix,'(a1,i0)') '.', iproc
# ifdef NETCDF_PARALLEL
       else
          WRITE (suffix,*) ''
       endif
# endif

       file_proc = trim(trim(file_proc)//adjustl(suffix))

# ifdef NETCDF_PARALLEL
       if(save_many) then
# endif
          istatus = nf90_create (file_proc, NF90_CLOBBER, ncid)
# ifdef NETCDF_PARALLEL
       else
          call barrier

          if(iproc .eq. 0) then
             open(unit=tmpunit, file=file_proc)
             close(unit=tmpunit, status='delete')
          end if

          call barrier
! If using netcdf version 4.1.2 or older replace NF90_MPIIO with NF90_CLOBBER
          istatus = nf90_create (file_proc, IOR(NF90_HDF5,NF90_MPIIO), ncid, comm=mp_comm, info=mp_info)
       end if
# endif

       if (istatus /= NF90_NOERR) then
          ierr = error_unit()
          write(ierr,*) "nf90_create error: ", nf90_strerror(istatus)
          goto 1
       end if

# ifdef NETCDF_PARALLEL
       if(.not.save_many) then
          istatus = nf90_put_att(ncid, NF90_GLOBAL, 'xyzs_layout', xyzs_layout)
          if (istatus /= NF90_NOERR) then
             ierr = error_unit()
             write(ierr,*) "nf90_put_attr error: ", nf90_strerror(istatus)
             goto 1
          end if
          istatus = nf90_put_att(ncid, NF90_GLOBAL, 'vms_layout', vms_layout)
          if (istatus /= NF90_NOERR) then
             ierr = error_unit()
             write(ierr,*) "nf90_put_attr error: ", nf90_strerror(istatus)
             goto 1
          end if
       endif
# endif

       if (n_elements > 0) then
          istatus = nf90_def_dim (ncid, "tube", ntubes, tubeid)
          if (istatus /= NF90_NOERR) then
             ierr = error_unit()
             write(ierr,*) "nf90_def_dim zed error: ", nf90_strerror(istatus)
             goto 1
          end if

          istatus = nf90_def_dim (ncid, "zed", 2*nzgrid+1, zedid)
          if (istatus /= NF90_NOERR) then
             ierr = error_unit()
             write(ierr,*) "nf90_def_dim zed error: ", nf90_strerror(istatus)
             goto 1
          end if

          istatus = nf90_def_dim (ncid, "vpa", nvpa, vpaid)
          if (istatus /= NF90_NOERR) then
             ierr = error_unit()
             write(ierr,*) "nf90_def_dim vpa error: ", nf90_strerror(istatus)
             goto 1
          end if

          istatus = nf90_def_dim (ncid, "mu", nmu, muid)
          if (istatus /= NF90_NOERR) then
             ierr = error_unit()
             write(ierr,*) "nf90_def_dim mu error: ", nf90_strerror(istatus)
             goto 1
          end if

# ifdef NETCDF_PARALLEL
          if(save_many) then
# endif
             istatus = nf90_def_dim (ncid, "glo", n_elements, gloid)
# ifdef NETCDF_PARALLEL
          else
             istatus = nf90_def_dim (ncid, "glo", total_elements, gloid)
          endif
# endif
          if (istatus /= NF90_NOERR) then
             ierr = error_unit()
             write(ierr,*) "nf90_def_dim glo error: ", nf90_strerror(istatus)
             goto 1
          end if

# ifdef NETCDF_PARALLEL
          if(save_many) then
# endif
             istatus = nf90_def_dim (ncid, "gvmulo", nvmulo_elements, gvmuloid)
# ifdef NETCDF_PARALLEL
          else
             istatus = nf90_def_dim (ncid, "gvmulo", total_vmulo_elements, gvmuloid)
          endif
# endif
          if (istatus /= NF90_NOERR) then
             ierr = error_unit()
             write(ierr,*) "nf90_def_dim gvmulo error: ", nf90_strerror(istatus)
             goto 1
          end if

          istatus = nf90_def_dim (ncid, "aky", naky, kyid)
          if (istatus /= NF90_NOERR) then
             ierr = error_unit()
             write(ierr,*) "nf90_def_dim aky error: ", nf90_strerror(istatus)
             goto 1
          end if

          istatus = nf90_def_dim (ncid, "akx", nakx, kxid)
          if (istatus /= NF90_NOERR) then
             ierr = error_unit()
             write(ierr,*) "nf90_def_dim akx error: ", nf90_strerror(istatus)
             goto 1
          end if
       end if

       if (netcdf_real == 0) netcdf_real = get_netcdf_code_precision()

       istatus = nf90_def_var (ncid, "t0", netcdf_real, t0id)
       if (istatus /= NF90_NOERR) then
          ierr = error_unit()
          write(ierr,*) "nf90_def_var t0 error: ", nf90_strerror(istatus)
          goto 1
       end if

       istatus = nf90_def_var (ncid, "istep0", nf90_int, istep0id)
       if (istatus /= NF90_NOERR) then
          ierr = error_unit()
          write(ierr,*) "nf90_def_var istep0 error: ", nf90_strerror(istatus)
          goto 1
       end if

       istatus = nf90_def_var (ncid, "delt0", netcdf_real, delt0id)
       if (istatus /= NF90_NOERR) then
          ierr = error_unit()
          write(ierr,*) "nf90_def_var delt0 error: ", nf90_strerror(istatus)
          goto 1
       end if

       if (n_elements > 0) then
          istatus = nf90_def_var (ncid, "gr", netcdf_real, &
               (/ vpaid, muid, gloid /), gr_id)
          if (istatus /= NF90_NOERR) then
             ierr = error_unit()
             write(ierr,*) "nf90_def_var g error: ", nf90_strerror(istatus)
             goto 1
          end if

          istatus = nf90_def_var (ncid, "gi", netcdf_real, &
               (/ vpaid, muid, gloid /), gi_id)
          if (istatus /= NF90_NOERR) then
             ierr = error_unit()
             write(ierr,*) "nf90_def_var g error: ", nf90_strerror(istatus)
             goto 1
          end if

          if (fphi > epsilon(0.)) then
             istatus = nf90_def_var (ncid, "phi_r", netcdf_real, &
                  (/ kyid, kxid, zedid, tubeid/), phir_id)
             if (istatus /= NF90_NOERR) then
                ierr = error_unit()
                write(ierr,*) "nf90_def_var phi error: ", nf90_strerror(istatus)
                goto 1
             end if

             istatus = nf90_def_var (ncid, "phi_i", netcdf_real, &
                  (/ kyid, kxid, zedid, tubeid /), phii_id)
             if (istatus /= NF90_NOERR) then
                ierr = error_unit()
                write(ierr,*) "nf90_def_var phi error: ", nf90_strerror(istatus)
                goto 1
             end if
          end if

          if (fapar > epsilon(0.)) then
             istatus = nf90_def_var (ncid, "apar_r", netcdf_real, &
                  (/ kyid, kxid, zedid, tubeid /), aparr_id)
             if (istatus /= NF90_NOERR) then
                ierr = error_unit()
                write(ierr,*) "nf90_def_var apar error: ", nf90_strerror(istatus)
                goto 1
             end if

             istatus = nf90_def_var (ncid, "apar_i", netcdf_real, &
                  (/ kyid, kxid, zedid, tubeid /), apari_id)
             if (istatus /= NF90_NOERR) then
                ierr = error_unit()
                write(ierr,*) "nf90_def_var apar error: ", nf90_strerror(istatus)
                goto 1
             end if
          end if

          write(*,*) "Maybe write bpar"
          ! Bob: write bpar info to netcdf
          if (fbpar > epsilon(0.)) then
            write(*,*) "try to write bpar"
             istatus = nf90_def_var (ncid, "bpar_r", netcdf_real, &
                  (/ kyid, kxid, zedid, tubeid /), bparr_id)
             if (istatus /= NF90_NOERR) then
                ierr = error_unit()
                write(ierr,*) "nf90_def_var bpar error: ", nf90_strerror(istatus)
                goto 1
             end if

             istatus = nf90_def_var (ncid, "bpar_i", netcdf_real, &
                  (/ kyid, kxid, zedid, tubeid /), bpari_id)
             if (istatus /= NF90_NOERR) then
                ierr = error_unit()
                write(ierr,*) "nf90_def_var bpar error: ", nf90_strerror(istatus)
                goto 1
             end if
          end if

          if (include_krook_operator) then
             istatus = nf90_def_var (ncid, "intkrook", netcdf_real, intkrook_id)
             if (istatus /= NF90_NOERR) then
                ierr = error_unit()
                write(ierr,*) "nf90_def_var intkrook error: ", nf90_strerror(istatus)
                goto 1
             end if

             istatus = nf90_def_var (ncid, "krookr", netcdf_real, &
                  (/ kxid, tubeid, gvmuloid /), krookr_id)
             if (istatus /= NF90_NOERR) then
                ierr = error_unit()
                write(ierr,*) "nf90_def_var apar error: ", nf90_strerror(istatus)
                goto 1
             end if

             istatus = nf90_def_var (ncid, "krooki", netcdf_real, &
                  (/ kxid, tubeid, gvmuloid /), krooki_id)
             if (istatus /= NF90_NOERR) then
                ierr = error_unit()
                write(ierr,*) "nf90_def_var krooki error: ", nf90_strerror(istatus)
                goto 1
             end if

          end if

          if (remove_zero_projection) then
             istatus = nf90_def_var (ncid, "intproj", netcdf_real, intproj_id)
             if (istatus /= NF90_NOERR) then
                ierr = error_unit()
                write(ierr,*) "nf90_def_var intproj error: ", nf90_strerror(istatus)
                goto 1
             end if

             istatus = nf90_def_var (ncid, "projr", netcdf_real, &
                  (/ tubeid, gvmuloid /), projr_id)
             if (istatus /= NF90_NOERR) then
                ierr = error_unit()
                write(ierr,*) "nf90_def_var projr error: ", nf90_strerror(istatus)
                goto 1
             end if

             istatus = nf90_def_var (ncid, "proji", netcdf_real, &
                  (/ tubeid, gvmuloid /), proji_id)
             if (istatus /= NF90_NOERR) then
                ierr = error_unit()
                write(ierr,*) "nf90_def_var proji error: ", nf90_strerror(istatus)
                goto 1
             end if

          end if

          if (prp_shear_enabled) then
             istatus = nf90_def_var (ncid, "shiftstate", netcdf_real,&
                (/ kyid /), shift_id)
             if (istatus /= NF90_NOERR) then
                ierr = error_unit()
                write(ierr,*) "nf90_def_var shiftstate error: ", nf90_strerror(istatus)
                goto 1
             end if
          endif

!           if (fbpar > epsilon(0.)) then
!              istatus = nf90_def_var (ncid, "bpar_r", netcdf_real, &
!                   (/ zedid, kxid, kyid /), bparr_id)
!              if (istatus /= NF90_NOERR) then
!                 ierr = error_unit()
!                 write(ierr,*) "nf90_def_var bparr error: ", nf90_strerror(istatus)
!                 goto 1
!              end if

!              istatus = nf90_def_var (ncid, "bpar_i", netcdf_real, &
!                   (/ zedid, kxid, kyid /), bpari_id)
!              if (istatus /= NF90_NOERR) then
!                 ierr = error_unit()
!                 write(ierr,*) "nf90_def_var bpari error: ", nf90_strerror(istatus)
!                 goto 1
!              end if
!           end if

       end if

! remove allocated conditional because we want to be able to restart
! using exb shear from a case which does not have exb shear (i.e.
! we need kx_shift variable defined in netcdf file even if no exb
! shear present in simulation) -- MAB + CMR
!       if (allocated(kx_shift)) then   ! MR begin
!       istatus = nf90_def_var (ncid, "kx_shift", netcdf_real, &
!            (/ kyid /), kx_shift_id)
!       if (istatus /= NF90_NOERR) then
!          ierr = error_unit()
!          write(ierr,*) "nf90_def_var kx_shift error: ", nf90_strerror(istatus)
!          goto 1
!       endif
!       endif   ! MR end

!    if (proc0) then
!      write (*,*) "Finished definitions"
    !      write (*,*) "List restart files"
    !      call system("echo 'defs' >> filelist.txt; ls nc/* >> filelist.txt;  ")
!    end if

       istatus = nf90_enddef (ncid)
       if (istatus /= NF90_NOERR) then
          ierr = error_unit()
          write (ierr,*) "nf90_enddef error: ", nf90_strerror(istatus)
          goto 1
       end if
    end if


    !!!-----------------------!!!
    !!!-----------------------!!!
    !!!-----------------------!!!

# ifdef NETCDF_PARALLEL
    if(save_many .or. iproc == 0) then
# endif

       istatus = nf90_put_var (ncid, delt0id, delt0)
       if (istatus /= NF90_NOERR) then
          ierr = error_unit()
          write (ierr,*) "nf90_put_var delt0 error: ", nf90_strerror(istatus)
          goto 1
       end if

       istatus = nf90_put_var (ncid, t0id, t0)
       if (istatus /= NF90_NOERR) then
          ierr = error_unit()
          write (ierr,*) "nf90_put_var t0 error: ", nf90_strerror(istatus)
          goto 1
       end if

       istatus = nf90_put_var (ncid, istep0id, istep0)
       if (istatus /= NF90_NOERR) then
          ierr = error_unit()
          write (ierr,*) "nf90_put_var istep0 error: ", nf90_strerror(istatus)
          goto 1
       end if

# ifdef NETCDF_PARALLEL
    endif
# endif

1   continue

    if (istatus /= NF90_NOERR) then
       i = nf90_close (ncid)
       return
    end if

    if (n_elements > 0) then

       if (.not. allocated(tmpr)) &
            allocate (tmpr(nvpa,nmu,kxkyz_lo%llim_proc:kxkyz_lo%ulim_alloc))

       tmpr = real(g)

# ifdef NETCDF_PARALLEL
       if(save_many) then
# endif
          istatus = nf90_put_var (ncid, gr_id, tmpr)
#ifdef NETCDF_PARALLEL
       else
          istatus = nf90_var_par_access(ncid, gr_id, NF90_COLLECTIVE)
          istatus = nf90_var_par_access(ncid, gi_id, NF90_COLLECTIVE)

          start_pos = (/1,1,kxkyz_lo%llim_proc+1/)
          counts = (/nvpa, nmu, n_elements/)

          istatus = nf90_put_var (ncid, gr_id, tmpr, start=start_pos, count=counts)
       endif
# endif

       if (istatus /= NF90_NOERR) call netcdf_error (istatus, ncid, gr_id)

       tmpr = aimag(g)
# ifdef NETCDF_PARALLEL
       if(save_many) then
# endif
          istatus = nf90_put_var (ncid, gi_id, tmpr)
#ifdef NETCDF_PARALLEL
       else
          istatus = nf90_put_var (ncid, gi_id, tmpr, start=start_pos, count=counts)
       endif
# endif

       if (istatus /= NF90_NOERR) call netcdf_error (istatus, ncid, gi_id)

# ifdef NETCDF_PARALLEL
       if(save_many .or. iproc == 0) then
# endif

          if (.not. allocated(ftmpr)) allocate (ftmpr(naky,nakx,2*nzgrid+1,ntubes))
          if (.not. allocated(ftmpi)) allocate (ftmpi(naky,nakx,2*nzgrid+1,ntubes))

          if (fphi > epsilon(0.)) then
             ftmpr = real(phi)
             istatus = nf90_put_var (ncid, phir_id, ftmpr)
             if (istatus /= NF90_NOERR) call netcdf_error (istatus, ncid, phir_id)

             ftmpi = aimag(phi)
             istatus = nf90_put_var (ncid, phii_id, ftmpi)
             if (istatus /= NF90_NOERR) call netcdf_error (istatus, ncid, phii_id)
          end if

          if (fapar > epsilon(0.)) then
             ftmpr = real(apar)
             istatus = nf90_put_var (ncid, aparr_id, ftmpr)
             if (istatus /= NF90_NOERR) call netcdf_error (istatus, ncid, aparr_id)

             ftmpi = aimag(apar)
             istatus = nf90_put_var (ncid, apari_id, ftmpi)
             if (istatus /= NF90_NOERR) call netcdf_error (istatus, ncid, apari_id)
          end if

          !!! Bob: Consider bpar

!           if (fbpar > epsilon(0.)) then
!              ftmpr = real(bparnew)
!              istatus = nf90_put_var (ncid, bparr_id, ftmpr)
!              if (istatus /= NF90_NOERR) call netcdf_error (istatus, ncid, bparr_id)

!              ftmpi = aimag(bparnew)
!              istatus = nf90_put_var (ncid, bpari_id, ftmpi)
!              if (istatus /= NF90_NOERR) call netcdf_error (istatus, ncid, bpari_id)
!           end if

# ifdef NETCDF_PARALLEL
       end if
# endif

       if (include_krook_operator) then
         if (.not. allocated(ktmpr)) &
           allocate (ktmpr(nakx,ntubes,vmu_lo%llim_proc:vmu_lo%ulim_alloc))
         if (.not. allocated(ktmpi)) &
           allocate (ktmpi(nakx,ntubes,vmu_lo%llim_proc:vmu_lo%ulim_alloc))

# ifdef NETCDF_PARALLEL
         if(save_many .or. iproc == 0) then
# endif

           istatus = nf90_put_var (ncid, intkrook_id, int_krook)
           if (istatus /= NF90_NOERR) then
             ierr = error_unit()
             write (ierr,*) "nf90_put_var int_krook error: ", nf90_strerror(istatus)
             goto 1
           end if

# ifdef NETCDF_PARALLEL
         endif
# endif

         ktmpr = real(g_krook)
         ktmpi = aimag(g_krook)

# ifdef NETCDF_PARALLEL
         if(save_many) then
# endif
           istatus = nf90_put_var (ncid, krookr_id, ktmpr)
#ifdef NETCDF_PARALLEL
         else
           istatus = nf90_var_par_access(ncid, krookr_id, NF90_COLLECTIVE)
           istatus = nf90_var_par_access(ncid, krooki_id, NF90_COLLECTIVE)

           start_pos = (/1,1,vmu_lo%llim_proc+1/)
           counts = (/nakx, ntubes, nvmulo_elements/)

           istatus = nf90_put_var (ncid, krookr_id, ktmpr, start=start_pos, count=counts)
         endif
# endif
         if (istatus /= NF90_NOERR) call netcdf_error (istatus, ncid, krookr_id)



# ifdef NETCDF_PARALLEL
         if(save_many) then
# endif
           istatus = nf90_put_var (ncid, krooki_id, ktmpi)
#ifdef NETCDF_PARALLEL
         else
           istatus = nf90_put_var (ncid, krooki_id, ktmpi, start=start_pos, count=counts)
         endif
# endif
         if (istatus /= NF90_NOERR) call netcdf_error (istatus, ncid, krooki_id)

       end if

       if (remove_zero_projection) then
         if (.not. allocated(ptmpr)) &
           allocate (ptmpr(nakx,ntubes,vmu_lo%llim_proc:vmu_lo%ulim_alloc))
         if (.not. allocated(ptmpi)) &
           allocate (ptmpi(nakx,ntubes,vmu_lo%llim_proc:vmu_lo%ulim_alloc))

# ifdef NETCDF_PARALLEL
         if(save_many .or. iproc == 0) then
# endif

           istatus = nf90_put_var (ncid, intproj_id, int_proj)
           if (istatus /= NF90_NOERR) then
             ierr = error_unit()
             write (ierr,*) "nf90_put_var int_proj error: ", nf90_strerror(istatus)
             goto 1
           end if

# ifdef NETCDF_PARALLEL
         endif
# endif

         ptmpr = real(g_proj)
         ptmpi = aimag(g_proj)

# ifdef NETCDF_PARALLEL
         if(save_many) then
# endif
           istatus = nf90_put_var (ncid, projr_id, ptmpr)
#ifdef NETCDF_PARALLEL
         else
           istatus = nf90_var_par_access(ncid, projr_id, NF90_COLLECTIVE)
           istatus = nf90_var_par_access(ncid, proji_id, NF90_COLLECTIVE)

           start_pos = (/1,1,vmu_lo%llim_proc+1/)
           counts = (/nakx,ntubes, nvmulo_elements/)

           istatus = nf90_put_var (ncid, projr_id, ptmpr, start=start_pos, count=counts)
         endif
# endif
         if (istatus /= NF90_NOERR) call netcdf_error (istatus, ncid, projr_id)

# ifdef NETCDF_PARALLEL
         if(save_many) then
# endif
           istatus = nf90_put_var (ncid, proji_id, ptmpi)
#ifdef NETCDF_PARALLEL
         else
           istatus = nf90_put_var (ncid, proji_id, ptmpi, start=start_pos, count=counts)
         endif
# endif
         if (istatus /= NF90_NOERR) call netcdf_error (istatus, ncid, proji_id)

       end if

       if (prp_shear_enabled) then
         istatus = nf90_put_var (ncid, shift_id, shift_state)
         if (istatus /= NF90_NOERR) call netcdf_error (istatus, ncid, shift_id)
       end if
    end if

    if (exit) then
       i = nf90_close (ncid)
       if (i /= NF90_NOERR) &
            call netcdf_error (istatus, message='nf90_close error')
    else
       i = nf90_sync (ncid)
       if (i /= NF90_NOERR) &
            call netcdf_error (istatus, message='nf90_sync error')
    end if

# else

    if (proc0) write (error_unit(),*) &
         'WARNING: stella_save_for_restart is called without netcdf library'

# endif

    if (allocated(tmpr))  deallocate (tmpr)
    if (allocated(tmpi))  deallocate (tmpi)
    if (allocated(ftmpr)) deallocate (ftmpr)
    if (allocated(ftmpi)) deallocate (ftmpi)
    if (allocated(ptmpr)) deallocate (ptmpr)
    if (allocated(ptmpi)) deallocate (ptmpi)
    if (allocated(ktmpr)) deallocate (ktmpr)
    if (allocated(ktmpi)) deallocate (ktmpi)

  end subroutine stella_save_for_restart



!!----------------------------------------------------------------------!!
!!----------------------------------------------------------------------!!
!!---Restart------------------------------------------------------------!!
!!----------------------------------------------------------------------!!
!!----------------------------------------------------------------------!!

  subroutine stella_restore_many (g, scale, istatus, fphi, fapar) ! Bob: Consider bpar
# ifdef NETCDF
    use mp, only: iproc
    use fields_arrays, only: phi, apar, shift_state
    use dist_fn_arrays, only: g_krook, g_proj
    use kt_grids, only: naky, nakx
# endif
    use zgrid, only: nzgrid, ntubes
    use vpamu_grids, only: nvpa, nmu
    use stella_layouts, only: kxkyz_lo, vmu_lo
    use file_utils, only: error_unit
    use species, only: nspec
    use dissipation, only: include_krook_operator, int_krook
    use dissipation, only: remove_zero_projection, int_proj
    use physics_flags, only: prp_shear_enabled

    implicit none

    complex, dimension (:,:,kxkyz_lo%llim_proc:), intent (out) :: g
    real, intent (in) :: scale
    integer, intent (out) :: istatus
    real, intent (in) :: fphi, fapar
# ifdef NETCDF
# ifdef NETCDF_PARALLEL
    integer, dimension(3) :: counts, start_pos
# endif
    character (306) :: file_proc
    character (10) :: suffix
    integer :: i, n_elements, nvmulo_elements, ierr
    real :: fac

    n_elements = kxkyz_lo%ulim_proc-kxkyz_lo%llim_proc+1
    nvmulo_elements = vmu_lo%ulim_proc-vmu_lo%llim_proc+1
    if (n_elements <= 0) return

    if (.not.initialized) then
!       initialized = .true.
       file_proc = trim(restart_file)

# ifdef NETCDF_PARALLEL
       if(read_many) then
# endif
          write (suffix,'(a1,i0)') '.', iproc
          file_proc = trim(trim(file_proc)//adjustl(suffix))
          istatus = nf90_open (file_proc, NF90_NOWRITE, ncid)
# ifdef NETCDF_PARALLEL
       else
! If using netcdf version 4.1.2 deleted NF90_MPIIO and the associated IOR
          istatus = nf90_open (file_proc, IOR(NF90_NOWRITE, NF90_MPIIO), ncid, comm=mp_comm, info=mp_info)
       endif
# endif

       if (istatus /= NF90_NOERR) call netcdf_error (istatus, file=file_proc)

       ! check precision
       if (netcdf_real == 0) netcdf_real = get_netcdf_code_precision()
       call check_netcdf_file_precision (ncid)

       istatus = nf90_inq_dimid (ncid, "tube", tubeid)
       if (istatus /= NF90_NOERR) call netcdf_error (istatus, dim='tube')

       istatus = nf90_inq_dimid (ncid, "zed", zedid)
       if (istatus /= NF90_NOERR) call netcdf_error (istatus, dim='zed')

       istatus = nf90_inq_dimid (ncid, "aky", kyid)
       if (istatus /= NF90_NOERR) call netcdf_error (istatus, dim='aky')

       istatus = nf90_inq_dimid (ncid, "akx", kxid)
       if (istatus /= NF90_NOERR) call netcdf_error (istatus, dim='akx')

       istatus = nf90_inq_dimid (ncid, "glo", gloid)
       if (istatus /= NF90_NOERR) call netcdf_error (istatus, dim='glo')

       istatus = nf90_inq_dimid (ncid, "gvmulo", gvmuloid)
       if (istatus /= NF90_NOERR) call netcdf_error (istatus, dim='gvmulo')

       istatus = nf90_inquire_dimension (ncid, tubeid, len=i)
       if (istatus /= NF90_NOERR) call netcdf_error (istatus, ncid, dimid=tubeid)
       if (i /= ntubes) write(*,*) 'Restart error: ntubes=? ',i,' : ',ntubes,' : ',iproc

       istatus = nf90_inquire_dimension (ncid, zedid, len=i)
       if (istatus /= NF90_NOERR) call netcdf_error (istatus, ncid, dimid=zedid)
       if (i /= 2*nzgrid + 1) write(*,*) 'Restart error: nzgrid=? ',i,' : ',nzgrid,' : ',iproc

       istatus = nf90_inquire_dimension (ncid, kyid, len=i)
       if (istatus /= NF90_NOERR) call netcdf_error (istatus, ncid, dimid=kyid)
       if (i /= naky) write(*,*) 'Restart error: naky=? ',i,' : ',naky,' : ',iproc

       istatus = nf90_inquire_dimension (ncid, kxid, len=i)
       if (istatus /= NF90_NOERR) call netcdf_error (istatus, ncid, dimid=kxid)
       if (i /= nakx) write(*,*) 'Restart error: nakx=? ',i,' : ',nakx,' : ',iproc

       istatus = nf90_inquire_dimension (ncid, gloid, len=i)
       if (istatus /= NF90_NOERR) call netcdf_error (istatus, ncid, dimid=gloid)
#ifdef NETCDF_PARALLEL
       if(read_many) then
#endif
          if (i /= kxkyz_lo%ulim_proc-kxkyz_lo%llim_proc+1) write(*,*) 'Restart error: glo=? ',i,' : ',iproc
#ifdef NETCDF_PARALLEL
       else
          if (i /= kxkyz_lo%ulim_world+1) write(*,*) 'Restart error: glo=? ',i,' : ',iproc
       endif
#endif
       istatus = nf90_inquire_dimension (ncid, gvmuloid, len=i)
       if (istatus /= NF90_NOERR) call netcdf_error (istatus, ncid, dimid=gvmuloid)
#ifdef NETCDF_PARALLEL
       if(read_many) then
#endif
          if (i /= vmu_lo%ulim_proc-vmu_lo%llim_proc+1) write(*,*) 'Restart error: gvmulo=? ',i,' : ',iproc
#ifdef NETCDF_PARALLEL
       else
          if (i /= vmu_lo%ulim_world+1) write(*,*) 'Restart error: gvmulo=? ',i,' : ',iproc
       endif
#endif

       if (fphi > epsilon(0.)) then
          istatus = nf90_inq_varid (ncid, "phi_r", phir_id)
          if (istatus /= NF90_NOERR) call netcdf_error (istatus, var='phi_r')

          istatus = nf90_inq_varid (ncid, "phi_i", phii_id)
          if (istatus /= NF90_NOERR) call netcdf_error (istatus, var='phi_i')
       end if

       if (fapar > epsilon(0.)) then
          istatus = nf90_inq_varid (ncid, "apar_r", aparr_id)
          if (istatus /= NF90_NOERR) call netcdf_error (istatus, var='apar_r')

          istatus = nf90_inq_varid (ncid, "apar_i", apari_id)
          if (istatus /= NF90_NOERR) call netcdf_error (istatus, var='apar_i')
       end if

       if(include_krook_operator) then
          istatus = nf90_inq_varid (ncid, "intkrook", intkrook_id)
          if (istatus /= NF90_NOERR) call netcdf_error (istatus, var='intkrook')

          istatus = nf90_inq_varid (ncid, "krookr", krookr_id)
          if (istatus /= NF90_NOERR) call netcdf_error (istatus, var='krookr')

          istatus = nf90_inq_varid (ncid, "krooki", krooki_id)
          if (istatus /= NF90_NOERR) call netcdf_error (istatus, var='krooki')

       endif

       if(remove_zero_projection) then
          istatus = nf90_inq_varid (ncid, "intproj", intproj_id)
          if (istatus /= NF90_NOERR) call netcdf_error (istatus, var='intproj')

          istatus = nf90_inq_varid (ncid, "projr", projr_id)
          if (istatus /= NF90_NOERR) call netcdf_error (istatus, var='projr')

          istatus = nf90_inq_varid (ncid, "proji", proji_id)
          if (istatus /= NF90_NOERR) call netcdf_error (istatus, var='proji')
       endif

       if(prp_shear_enabled) then
          istatus = nf90_inq_varid (ncid, "shiftstate", shift_id)
          if (istatus /= NF90_NOERR) call netcdf_error (istatus, var='shiftstate')
       endif

!        if (fbpar > epsilon(0.)) then
!           istatus = nf90_inq_varid (ncid, "bpar_r", bparr_id)
!           if (istatus /= NF90_NOERR) call netcdf_error (istatus, var='bpar_r')

!           istatus = nf90_inq_varid (ncid, "bpar_i", bpari_id)
!           if (istatus /= NF90_NOERR) call netcdf_error (istatus, var='bpar_i')
!        end if

!       if (allocated(kx_shift)) then   ! MR begin
!          istatus = nf90_inq_varid (ncid, "kx_shift", kx_shift_id)
!          if (istatus /= NF90_NOERR) call netcdf_error (istatus, var='kx_shift')
!       endif   ! MR end

       istatus = nf90_inq_varid (ncid, "gr", gr_id)
       if (istatus /= NF90_NOERR) call netcdf_error (istatus, var='gr')

       istatus = nf90_inq_varid (ncid, "gi", gi_id)
       if (istatus /= NF90_NOERR) call netcdf_error (istatus, var='gi')
    end if

    if (.not. allocated(tmpr)) &
         allocate (tmpr(nvpa,nmu,kxkyz_lo%llim_proc:kxkyz_lo%ulim_alloc))
    if (.not. allocated(tmpi)) &
         allocate (tmpi(nvpa,nmu,kxkyz_lo%llim_proc:kxkyz_lo%ulim_alloc))

    tmpr = 0.; tmpi = 0.
# ifdef NETCDF_PARALLEL
    if(read_many) then
# endif
       istatus = nf90_get_var (ncid, gr_id, tmpr)
#ifdef NETCDF_PARALLEL
    else
       start_pos = (/1,1,kxkyz_lo%llim_proc+1/)
       counts = (/nvpa, nmu, n_elements/)
       istatus = nf90_get_var (ncid, gr_id, tmpr, start=start_pos, count=counts)
    end if
# endif

   if (istatus /= NF90_NOERR) call netcdf_error (istatus, ncid, gr_id)

# ifdef NETCDF_PARALLEL
    if(read_many) then
# endif
       istatus = nf90_get_var (ncid, gi_id, tmpi)
#ifdef NETCDF_PARALLEL
    else
       istatus = nf90_get_var (ncid, gi_id, tmpi, start=start_pos, count=counts)
    end if
# endif

    if (istatus /= NF90_NOERR) call netcdf_error (istatus, ncid, gi_id)

    g = cmplx(tmpr, tmpi)

    if (.not. allocated(ftmpr)) allocate (ftmpr(naky,nakx,2*nzgrid+1,ntubes))
    if (.not. allocated(ftmpi)) allocate (ftmpi(naky,nakx,2*nzgrid+1,ntubes))

!    if (allocated(kx_shift)) then   ! MR begin
!       if (.not. allocated(stmp)) allocate (stmp(naky))   ! MR
!       istatus = nf90_get_var (ncid, kx_shift_id, stmp)
!       if (istatus /= NF90_NOERR) call netcdf_error (istatus, ncid, kx_shift_id)
!       kx_shift = stmp
!    endif   ! MR end

    if (fphi > epsilon(0.)) then
       istatus = nf90_get_var (ncid, phir_id, ftmpr)
       if (istatus /= NF90_NOERR) call netcdf_error (istatus, ncid, phir_id)

       istatus = nf90_get_var (ncid, phii_id, ftmpi)
       if (istatus /= NF90_NOERR) call netcdf_error (istatus, ncid, phii_id)

       phi = cmplx(ftmpr, ftmpi)
    end if

    if (fapar > epsilon(0.)) then
       istatus = nf90_get_var (ncid, aparr_id, ftmpr)
       if (istatus /= NF90_NOERR) call netcdf_error (istatus, ncid, aparr_id)

       istatus = nf90_get_var (ncid, apari_id, ftmpi)
       if (istatus /= NF90_NOERR) call netcdf_error (istatus, ncid, apari_id)

       apar = cmplx(ftmpr, ftmpi)
    end if

    if(include_krook_operator) then
      if (.not. allocated(ktmpr)) &
        allocate (ktmpr(nakx,ntubes,vmu_lo%llim_proc:vmu_lo%ulim_alloc))
      if (.not. allocated(ktmpi)) &
        allocate (ktmpi(nakx,ntubes,vmu_lo%llim_proc:vmu_lo%ulim_alloc))

      istatus = nf90_get_var (ncid, intkrook_id, int_krook)
      if (istatus /= NF90_NOERR) call netcdf_error (istatus, ncid, intkrook_id)

      ktmpr = 0.; ktmpi = 0.
# ifdef NETCDF_PARALLEL
      if(read_many) then
# endif
        istatus = nf90_get_var (ncid, krookr_id, ktmpr)
#ifdef NETCDF_PARALLEL
      else
        start_pos = (/1,1,vmu_lo%llim_proc+1/)
        counts = (/nakx, ntubes, nvmulo_elements/)
        istatus = nf90_get_var (ncid, krookr_id, ktmpr, start=start_pos, count=counts)
      end if
# endif

       if (istatus /= NF90_NOERR) call netcdf_error (istatus, ncid, krookr_id)

# ifdef NETCDF_PARALLEL
      if(read_many) then
# endif
        istatus = nf90_get_var (ncid, krooki_id, ktmpi)
#ifdef NETCDF_PARALLEL
      else
        istatus = nf90_get_var (ncid, krooki_id, ktmpi, start=start_pos, count=counts)
      end if
# endif

      if (istatus /= NF90_NOERR) call netcdf_error (istatus, ncid, krooki_id)

      g_krook = cmplx(ktmpr, ktmpi)

    endif

    if(remove_zero_projection) then
      if (.not. allocated(ptmpr)) &
        allocate (ptmpr(nakx,ntubes,vmu_lo%llim_proc:vmu_lo%ulim_alloc))
      if (.not. allocated(ptmpi)) &
        allocate (ptmpi(nakx,ntubes,vmu_lo%llim_proc:vmu_lo%ulim_alloc))

      istatus = nf90_get_var (ncid, intproj_id, int_proj)
      if (istatus /= NF90_NOERR) call netcdf_error (istatus, ncid, intproj_id)

      ptmpr = 0.; ptmpi = 0.
# ifdef NETCDF_PARALLEL
      if(read_many) then
# endif
        istatus = nf90_get_var (ncid, projr_id, ptmpr)
#ifdef NETCDF_PARALLEL
      else
        start_pos = (/1,1,vmu_lo%llim_proc+1/)
        counts = (/nakx,ntubes, nvmulo_elements/)
        istatus = nf90_get_var (ncid, projr_id, ptmpr, start=start_pos, count=counts)
      end if
# endif

       if (istatus /= NF90_NOERR) call netcdf_error (istatus, ncid, projr_id)

# ifdef NETCDF_PARALLEL
      if(read_many) then
# endif
        istatus = nf90_get_var (ncid, proji_id, ptmpi)
#ifdef NETCDF_PARALLEL
      else
        istatus = nf90_get_var (ncid, proji_id, ptmpi, start=start_pos, count=counts)
      end if
# endif

      if (istatus /= NF90_NOERR) call netcdf_error (istatus, ncid, proji_id)

      g_proj = cmplx(ptmpr, ptmpi)

    endif

    if(prp_shear_enabled) then
      istatus = nf90_get_var (ncid, shift_id, shift_state)
      if (istatus /= NF90_NOERR) call netcdf_error (istatus, ncid, shift_id)
    endif


!     if (fbpar > epsilon(0.)) then
!        istatus = nf90_get_var (ncid, bparr_id, ftmpr)
!        if (istatus /= NF90_NOERR) call netcdf_error (istatus, ncid, bparr_id)

!        istatus = nf90_get_var (ncid, bpari_id, ftmpi)
!        if (istatus /= NF90_NOERR) call netcdf_error (istatus, ncid, bpari_id)

!        bparold = 0.
!        bparnew = cmplx(ftmpr, ftmpi)
!     end if

    if (scale > 0.) then
       g = g*scale
       phi = phi*scale
       apar = apar*scale
       if(include_krook_operator) g_krook = g_krook*scale
       if(remove_zero_projection) g_proj = g_proj*scale
    else
       fac = - scale/(maxval(abs(phi)))
       g = g*fac
       phi = phi*fac
       apar = apar*fac
       if(include_krook_operator) g_krook = g_krook*fac
       if(remove_zero_projection) g_proj = g_proj*fac
    end if

    ! RN 2008/05/23: this was commented out. why? HJL 2013/05/15 Because it stops future writing to the file
!    istatus = nf90_close (ncid)
    if (istatus /= NF90_NOERR) then
       ierr = error_unit()
       write(ierr,*) "nf90_close error: ", nf90_strerror(istatus),' ',iproc
    end if

# else

    write (error_unit(),*) &
         'ERROR: stella_restore_many is called without netcdf'

# endif

    if (allocated(tmpr))  deallocate (tmpr)
    if (allocated(tmpi))  deallocate (tmpi)
    if (allocated(ftmpr)) deallocate (ftmpr)
    if (allocated(ftmpi)) deallocate (ftmpi)
    if (allocated(ptmpr)) deallocate (ptmpr)
    if (allocated(ptmpi)) deallocate (ptmpi)
    if (allocated(ktmpr)) deallocate (ktmpr)
    if (allocated(ktmpi)) deallocate (ktmpi)

  end subroutine stella_restore_many


  subroutine init_save (file)

    character(300), intent (in) :: file

    restart_file = file

  end subroutine init_save

  subroutine init_dt (delt0, istatus)

# ifdef NETCDF
    use mp, only: proc0, broadcast
    use file_utils, only: error_unit
# endif
    implicit none
    real, intent (in out) :: delt0
    integer, intent (out) :: istatus
# ifdef NETCDF
    character (306) :: file_proc

    if (proc0) then

       if (.not. initialized) then

# ifdef NETCDF_PARALLEL
          if(read_many) then
# endif
             file_proc=trim(trim(restart_file)//'.0')
# ifdef NETCDF_PARALLEL
          else
             file_proc=trim(trim(restart_file))
          end if
# endif

          istatus = nf90_open (file_proc, NF90_NOWRITE, ncid)
          if (istatus /= NF90_NOERR) call netcdf_error (istatus,file=file_proc)

          istatus = nf90_inq_varid (ncid, "delt0", delt0id)
          if (istatus /= NF90_NOERR) call netcdf_error (istatus, var='delt0')
       end if

       istatus = nf90_get_var (ncid, delt0id, delt0)

       if (istatus /= NF90_NOERR) then
          call netcdf_error (istatus, ncid, delt0id, message=' in init_dt')
          delt0 = -1.
       endif

       if (.not.initialized) istatus = nf90_close (ncid)
    endif

    call broadcast (istatus)
    call broadcast (delt0)

# endif

  end subroutine init_dt

  subroutine init_tstart (tstart, istep0, istatus)

# ifdef NETCDF
    use mp, only: proc0, broadcast
    use file_utils, only: error_unit
# endif
    implicit none
    real, intent (in out) :: tstart
    integer, intent (out) :: istep0
    integer, intent (out) :: istatus
# ifdef NETCDF
    character (306) :: file_proc

    if (proc0) then
# ifdef NETCDF_PARALLEL
       if(read_many) then
# endif
          file_proc=trim(trim(restart_file)//'.0')
# ifdef NETCDF_PARALLEL
       else
          file_proc=trim(trim(restart_file))
       end if
# endif

       if (.not.initialized) then

          istatus = nf90_open (file_proc, NF90_NOWRITE, ncid)
          if (istatus /= NF90_NOERR) call netcdf_error (istatus, file=file_proc)
       end if

       istatus = nf90_inq_varid (ncid, "t0", t0id)
       if (istatus /= NF90_NOERR) call netcdf_error (istatus, var='t0')

       istatus = nf90_get_var (ncid, t0id, tstart)
       if (istatus /= NF90_NOERR) then
          call netcdf_error (istatus, ncid, t0id, message=' in init_tstart')
          tstart = -1.
       end if

       istatus = nf90_inq_varid (ncid, "istep0", istep0id)
       if (istatus /= NF90_NOERR) call netcdf_error (istatus, var='istep0')

       istatus = nf90_get_var (ncid, istep0id, istep0)
       if (istatus /= NF90_NOERR) then
          call netcdf_error (istatus, ncid, istep0id, message=' in init_tstart')
          istep0 = -1
       end if

       if (.not.initialized) istatus = nf90_close (ncid)

    endif

    call broadcast (istatus)
    call broadcast (istep0)
    call broadcast (tstart)

# endif

  end subroutine init_tstart

  subroutine finish_save

    if (allocated(tmpr))  deallocate (tmpr)
    if (allocated(tmpi))  deallocate (tmpi)
    if (allocated(ftmpr)) deallocate (ftmpr)
    if (allocated(ftmpi)) deallocate (ftmpi)

  end subroutine finish_save

end module stella_save
