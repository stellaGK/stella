!###############################################################################
!################################ GVEC GEOMETRY ################################
!###############################################################################
!
! Reads geometry arrays pre-computed by the Python script gvec_to_stella.py
! from a netCDF file and provides them to stella's geometry module.
!
! The file was produced as follows:
!   python gvec_to_stella.py  STATE_FILE  --s 0.5  --nalpha 1  --nzgrid 128
!
! All normalisation conventions match stella's VMEC reader.  See
! geometry_vmec_read_netCDF_file.f90 for definitions of the quantities.
!
!###############################################################################
module geometry_gvec_read_netCDF_file

   use debug_flags, only: debug => geometry_debug

   implicit none

   public :: read_gvec_geometry

   private

contains

   !============================================================================
   !================== READ GVEC GEOMETRY NETCDF FILE ==========================
   !============================================================================
   ! Reads the gvec_stella_geometry*.nc file produced by gvec_to_stella.py
   ! and returns all geometry arrays needed by stella.
   !
   ! The output arrays use the same normalisation as calculate_vmec_geometry.
   !============================================================================
   subroutine read_gvec_geometry( &
      ! Input
      filename, nalpha, nzgrid_file, verbose, &
      ! Scalar outputs
      s, safety_factor_q, shat, L_reference, B_reference, nfp_out, &
      sign_toroidal_flux, nfp_sim_out, dzdzeta, alpha0_file, zeta_center_file, &
      ! 1-D outputs
      alpha, zeta, &
      ! 2-D outputs  (nalpha, -nzgrid_file:nzgrid_file)
      theta, grho, bmag, b_dot_gradzeta, &
      B_sub_theta, B_sub_zeta, &
      grad_alpha_grad_alpha, grad_alpha_grad_psit, grad_psit_grad_psit, &
      gds23_psitalpha, gds24_psitalpha, gds25_psitalpha, gds26_psitalpha, &
      gbdrift_alpha, B_times_gradB_dot_gradpsit, &
      gbdrift_kappa_alpha, B_times_kappa_dot_gradpsit, &
      psit_displacement_fac, &
      gradzeta_gradalpha_R2overB2, gradzeta_gradpsit_R2overB2, &
      b_dot_gradzeta_RR, &
      ierr)

#ifdef NETCDF
      use netcdf, only: nf90_open, nf90_close, nf90_nowrite, nf90_noerr
      use netcdf, only: nf90_inq_varid, nf90_get_var, nf90_inq_dimid
      use netcdf, only: nf90_inquire_dimension, nf90_inquire_variable
#endif
      use mp, only: mp_abort
      use debug_flags, only: print_extra_info_to_terminal

      implicit none

      ! Input
      character(*), intent(in) :: filename
      integer,      intent(in) :: nalpha, nzgrid_file
      logical,      intent(in) :: verbose

      ! Scalar outputs
      integer, intent(out) :: sign_toroidal_flux, nfp_out, ierr
      integer, intent(out) :: nfp_sim_out     ! number of simulated field periods
      real,    intent(out) :: s, safety_factor_q, shat
      real,    intent(out) :: L_reference, B_reference, dzdzeta
      real,    intent(out) :: alpha0_file, zeta_center_file  ! for namelist validation

      ! 1-D outputs
      real, dimension(:),            intent(out) :: alpha
      real, dimension(-nzgrid_file:), intent(out) :: zeta

      ! 2-D outputs
      real, dimension(:, -nzgrid_file:), intent(out) :: theta, grho, bmag
      real, dimension(:, -nzgrid_file:), intent(out) :: b_dot_gradzeta
      real, dimension(:, -nzgrid_file:), intent(out) :: B_sub_theta, B_sub_zeta
      real, dimension(:, -nzgrid_file:), intent(out) :: grad_alpha_grad_alpha
      real, dimension(:, -nzgrid_file:), intent(out) :: grad_alpha_grad_psit
      real, dimension(:, -nzgrid_file:), intent(out) :: grad_psit_grad_psit
      real, dimension(:, -nzgrid_file:), intent(out) :: gds23_psitalpha
      real, dimension(:, -nzgrid_file:), intent(out) :: gds24_psitalpha
      real, dimension(:, -nzgrid_file:), intent(out) :: gds25_psitalpha
      real, dimension(:, -nzgrid_file:), intent(out) :: gds26_psitalpha
      real, dimension(:, -nzgrid_file:), intent(out) :: gbdrift_alpha
      real, dimension(:, -nzgrid_file:), intent(out) :: B_times_gradB_dot_gradpsit
      real, dimension(:, -nzgrid_file:), intent(out) :: gbdrift_kappa_alpha
      real, dimension(:, -nzgrid_file:), intent(out) :: B_times_kappa_dot_gradpsit
      real, dimension(:, -nzgrid_file:), intent(out) :: psit_displacement_fac
      real, dimension(:, -nzgrid_file:), intent(out) :: gradzeta_gradalpha_R2overB2
      real, dimension(:, -nzgrid_file:), intent(out) :: gradzeta_gradpsit_R2overB2
      real, dimension(:, -nzgrid_file:), intent(out) :: b_dot_gradzeta_RR

#ifdef NETCDF
      ! Local variables
      integer :: ncid, varid, dimid
      integer :: nz_file_check, nalpha_check, nzgrid_check
      integer :: nfp_tmp, sign_tmp
      ! Temporary flat 2-D buffer (C-order from netCDF): (nalpha, nz_file)
      real, allocatable :: buf2d(:, :)   ! shape (nalpha, 2*nzgrid_file+1)
      real(kind=8), allocatable :: buf1d(:)
      integer :: nz_file

      ierr = 0
      nz_file = 2 * nzgrid_file + 1

      if (debug) write (*, *) 'geometry_gvec_read_netCDF_file::read_gvec_geometry'

      if (verbose .and. print_extra_info_to_terminal) then
         write (*, '(A)') " "
         write (*, '(A)') "######################################################"
         write (*, '(A)') "               GVEC GEOMETRY"
         write (*, '(A)') "######################################################"
         write (*, '(A)') " "
         write (*, '(2A)') "  Reading GVEC geometry file: ", trim(filename)
      end if

      ! Open file
      call check(nf90_open(trim(filename), nf90_nowrite, ncid), "open", filename, ierr)
      if (ierr /= 0) then
         write (*, '(2A)') "  ERROR: cannot open GVEC geometry file: ", trim(filename)
         call mp_abort("read_gvec_geometry: file not found")
      end if

      ! Validate dimensions
      call check(nf90_inq_dimid(ncid, "nalpha", dimid), "dimid nalpha", filename, ierr)
      call check(nf90_inquire_dimension(ncid, dimid, len=nalpha_check), "dim nalpha", filename, ierr)
      if (nalpha_check /= nalpha) then
         write (*, '(A,2I6)') "  ERROR: nalpha in file vs stella: ", nalpha_check, nalpha
         call mp_abort("read_gvec_geometry: nalpha mismatch - re-run gvec_to_stella.py with --nalpha "//to_str(nalpha))
      end if

      call check(nf90_inq_dimid(ncid, "nz_file", dimid), "dimid nz_file", filename, ierr)
      call check(nf90_inquire_dimension(ncid, dimid, len=nz_file_check), "dim nz_file", filename, ierr)
      if (nz_file_check /= nz_file) then
         write (*, '(A,2I6)') "  ERROR: nz_file in file vs expected: ", nz_file_check, nz_file
         call mp_abort("read_gvec_geometry: nzgrid_file mismatch - re-run gvec_to_stella.py with --nzgrid_file "//to_str(nzgrid_file))
      end if

      ! Scalar integers
      call read_int_scalar(ncid, "nfp",          nfp_tmp,  filename, ierr); nfp_out = nfp_tmp
      call read_int_scalar(ncid, "sign_torflux", sign_tmp, filename, ierr); sign_toroidal_flux = sign_tmp
      call read_int_scalar(ncid, "nfp_sim",      nfp_sim_out, filename, ierr)

      ! Scalar reals
      call read_dbl_scalar(ncid, "s",               s,                filename, ierr)
      call read_dbl_scalar(ncid, "safety_factor_q", safety_factor_q,  filename, ierr)
      call read_dbl_scalar(ncid, "shat",             shat,             filename, ierr)
      call read_dbl_scalar(ncid, "L_reference",      L_reference,      filename, ierr)
      call read_dbl_scalar(ncid, "B_reference",      B_reference,      filename, ierr)
      call read_dbl_scalar(ncid, "dzdzeta",          dzdzeta,          filename, ierr)
      call read_dbl_scalar(ncid, "alpha0",           alpha0_file,      filename, ierr)
      call read_dbl_scalar(ncid, "zeta_center",      zeta_center_file, filename, ierr)

      ! 1-D arrays
      allocate (buf1d(nalpha))
      call read_1d(ncid, "alpha", buf1d, filename, ierr)
      alpha(:) = real(buf1d(:))
      deallocate (buf1d)

      allocate (buf1d(nz_file))
      call read_1d(ncid, "zeta", buf1d, filename, ierr)
      zeta(-nzgrid_file:nzgrid_file) = real(buf1d(:))
      deallocate (buf1d)

      ! 2-D arrays
      ! netCDF stores as (nalpha, nz_file) in C order = Fortran (nz_file, nalpha)
      ! We read into a temporary (nalpha, nz_file) buffer then assign.
      allocate (buf2d(nalpha, nz_file))

      call read_2d_and_store(ncid, "theta",                      buf2d, theta,                      nzgrid_file, filename, ierr)
      call read_2d_and_store(ncid, "grho",                       buf2d, grho,                       nzgrid_file, filename, ierr)
      call read_2d_and_store(ncid, "bmag",                       buf2d, bmag,                       nzgrid_file, filename, ierr)
      call read_2d_and_store(ncid, "b_dot_gradzeta",             buf2d, b_dot_gradzeta,             nzgrid_file, filename, ierr)
      call read_2d_and_store(ncid, "B_sub_theta",                buf2d, B_sub_theta,                nzgrid_file, filename, ierr)
      call read_2d_and_store(ncid, "B_sub_zeta",                 buf2d, B_sub_zeta,                 nzgrid_file, filename, ierr)
      call read_2d_and_store(ncid, "grad_alpha_grad_alpha",      buf2d, grad_alpha_grad_alpha,      nzgrid_file, filename, ierr)
      call read_2d_and_store(ncid, "grad_alpha_grad_psit",       buf2d, grad_alpha_grad_psit,       nzgrid_file, filename, ierr)
      call read_2d_and_store(ncid, "grad_psit_grad_psit",        buf2d, grad_psit_grad_psit,        nzgrid_file, filename, ierr)
      call read_2d_and_store(ncid, "gds23",                      buf2d, gds23_psitalpha,            nzgrid_file, filename, ierr)
      call read_2d_and_store(ncid, "gds24",                      buf2d, gds24_psitalpha,            nzgrid_file, filename, ierr)
      call read_2d_and_store(ncid, "gds25",                      buf2d, gds25_psitalpha,            nzgrid_file, filename, ierr)
      call read_2d_and_store(ncid, "gds26",                      buf2d, gds26_psitalpha,            nzgrid_file, filename, ierr)
      call read_2d_and_store(ncid, "gbdrift_alpha",              buf2d, gbdrift_alpha,              nzgrid_file, filename, ierr)
      call read_2d_and_store(ncid, "B_times_gradB_dot_gradpsit", buf2d, B_times_gradB_dot_gradpsit, nzgrid_file, filename, ierr)
      call read_2d_and_store(ncid, "gbdrift_kappa_alpha",        buf2d, gbdrift_kappa_alpha,        nzgrid_file, filename, ierr)
      call read_2d_and_store(ncid, "B_times_kappa_dot_gradpsit", buf2d, B_times_kappa_dot_gradpsit, nzgrid_file, filename, ierr)
      call read_2d_and_store(ncid, "psit_displacement_fac",      buf2d, psit_displacement_fac,      nzgrid_file, filename, ierr)
      call read_2d_and_store(ncid, "gradzeta_gradalpha_R2overB2",buf2d, gradzeta_gradalpha_R2overB2,nzgrid_file, filename, ierr)
      call read_2d_and_store(ncid, "gradzeta_gradpsit_R2overB2", buf2d, gradzeta_gradpsit_R2overB2, nzgrid_file, filename, ierr)
      call read_2d_and_store(ncid, "b_dot_gradzeta_RR",          buf2d, b_dot_gradzeta_RR,          nzgrid_file, filename, ierr)

      deallocate (buf2d)

      call check(nf90_close(ncid), "close", filename, ierr)

      if (verbose .and. print_extra_info_to_terminal) then
         write (*, '(A,F8.5)') "    s             = ", s
         write (*, '(A,F8.5)') "    q             = ", safety_factor_q
         write (*, '(A,F8.5)') "    shat          = ", shat
         write (*, '(A,F8.5)') "    L_reference   = ", L_reference
         write (*, '(A,F8.5)') "    B_reference   = ", B_reference
         write (*, '(A,I3)' )  "    nfp           = ", nfp_out
         write (*, '(A,I3)' )  "    sign_torflux  = ", sign_toroidal_flux
         write (*, *) " "
      end if

#else
      ierr = -1
      call mp_abort("read_gvec_geometry: stella was compiled without netCDF support")
#endif

   end subroutine read_gvec_geometry


#ifdef NETCDF
   !============================================================================
   ! Internal helper: read a 2-D variable (nalpha, nz_file) into a
   ! Fortran array with index range (nalpha, -nzgrid_file:nzgrid_file).
   !============================================================================
   subroutine read_2d_and_store(ncid, varname, buf2d, out2d, nzgrid_file, filename, ierr)
      use netcdf, only: nf90_inq_varid, nf90_get_var, nf90_noerr
      implicit none
      integer,          intent(in)    :: ncid, nzgrid_file
      character(*),     intent(in)    :: varname, filename
      real,             intent(inout) :: buf2d(:, :)   ! (nalpha, nz_file) workspace
      real, dimension(:, -nzgrid_file:), intent(out) :: out2d
      integer,          intent(inout) :: ierr

      integer :: varid, stat, ia, iz, nalpha_loc, nz_loc
      double precision, allocatable :: dbl_buf(:, :)

      nalpha_loc = size(buf2d, 1)
      nz_loc     = size(buf2d, 2)   ! = 2*nzgrid_file + 1

      ! NOTE: netCDF stores the variable with dimensions (nalpha, nz_file) as
      ! declared by gvec_to_stella.py (C order). The Fortran netCDF API presents
      ! dimensions in REVERSE order, so the receiving buffer must be
      ! (nz_file, nalpha). Declaring it (nalpha, nz_file) makes every read fail
      ! with "Start+count exceeds dimension bound" - including for nalpha = 1.
      allocate (dbl_buf(nz_loc, nalpha_loc))

      ! Only the covariant components of B, the displacement factor and the factors of the
      ! momentum flux may be missing from the file (they are then set to zero); the other
      ! arrays define the geometry, so a missing one is an error
      stat = nf90_inq_varid(ncid, varname, varid)
      if (stat /= nf90_noerr) then
         select case (trim(varname))
         case ('B_sub_theta', 'B_sub_zeta', 'psit_displacement_fac', 'gradzeta_gradalpha_R2overB2', &
               'gradzeta_gradpsit_R2overB2', 'b_dot_gradzeta_RR')
            write (*, '(3A)') "  WARNING: variable '", trim(varname), "' not found in GVEC geometry file, it is set to zero"
         case default
            write (*, '(3A)') "  ERROR: variable '", trim(varname), "' not found in GVEC geometry file"
            ierr = ierr + 1
         end select
         out2d = 0.0
         deallocate(dbl_buf)
         return
      end if

      call check_stat(nf90_get_var(ncid, varid, dbl_buf), varname, filename, ierr)

      do ia = 1, nalpha_loc
         do iz = 1, nz_loc
            out2d(ia, iz - nzgrid_file - 1) = real(dbl_buf(iz, ia))
         end do
      end do

      deallocate (dbl_buf)
   end subroutine read_2d_and_store


   subroutine read_1d(ncid, varname, buf, filename, ierr)
      use netcdf, only: nf90_inq_varid, nf90_get_var
      implicit none
      integer,      intent(in)    :: ncid
      character(*), intent(in)    :: varname, filename
      real(kind=8), intent(out)   :: buf(:)
      integer,      intent(inout) :: ierr
      integer :: varid
      call check_stat(nf90_inq_varid(ncid, varname, varid),         varname, filename, ierr)
      call check_stat(nf90_get_var(ncid, varid, buf),               varname, filename, ierr)
   end subroutine read_1d


   subroutine read_int_scalar(ncid, varname, val, filename, ierr)
      use netcdf, only: nf90_inq_varid, nf90_get_var
      implicit none
      integer,      intent(in)    :: ncid
      character(*), intent(in)    :: varname, filename
      integer,      intent(out)   :: val
      integer,      intent(inout) :: ierr
      integer :: varid
      call check_stat(nf90_inq_varid(ncid, varname, varid), varname, filename, ierr)
      call check_stat(nf90_get_var(ncid, varid, val),       varname, filename, ierr)
   end subroutine read_int_scalar


   subroutine read_dbl_scalar(ncid, varname, val, filename, ierr)
      use netcdf, only: nf90_inq_varid, nf90_get_var
      implicit none
      integer,      intent(in)    :: ncid
      character(*), intent(in)    :: varname, filename
      real,         intent(out)   :: val
      integer,      intent(inout) :: ierr
      integer       :: varid
      real(kind=8)  :: tmp
      call check_stat(nf90_inq_varid(ncid, varname, varid), varname, filename, ierr)
      call check_stat(nf90_get_var(ncid, varid, tmp),       varname, filename, ierr)
      val = real(tmp)
   end subroutine read_dbl_scalar


   subroutine check(stat, context, filename, ierr)
      use netcdf, only: nf90_noerr, nf90_strerror
      use mp,     only: mp_abort
      implicit none
      integer,      intent(in)    :: stat
      character(*), intent(in)    :: context, filename
      integer,      intent(inout) :: ierr
      if (stat /= nf90_noerr) then
         write (*, '(4A)') "  netCDF error in '", trim(context), "' for file: ", trim(filename)
         write (*, '(A)')  nf90_strerror(stat)
         ierr = ierr + 1
      end if
   end subroutine check


   subroutine check_stat(stat, varname, filename, ierr)
      use netcdf, only: nf90_noerr, nf90_strerror
      implicit none
      integer,      intent(in)    :: stat
      character(*), intent(in)    :: varname, filename
      integer,      intent(inout) :: ierr
      if (stat /= nf90_noerr) then
         write (*, '(4A)') "  netCDF error reading '", trim(varname), "' from: ", trim(filename)
         write (*, '(A)')  nf90_strerror(stat)
         ierr = ierr + 1
      end if
   end subroutine check_stat


   ! Integer → string (simple, avoids format dependency)
   function to_str(i) result(s)
      implicit none
      integer, intent(in) :: i
      character(len=12) :: s
      write (s, '(I12)') i
      s = adjustl(s)
   end function to_str

#endif

end module geometry_gvec_read_netCDF_file
