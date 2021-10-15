MODULE ezcdf_opncls
!DEC$ IF DEFINED (NETCDF)
   INTERFACE cdfOpn
      MODULE PROCEDURE ezcdf_open
   END INTERFACE

   INTERFACE cdf_open
      MODULE PROCEDURE ezcdf_open
   END INTERFACE

   INTERFACE cdfCls
      MODULE PROCEDURE ezcdf_close
   END INTERFACE

   INTERFACE cdf_close
      MODULE PROCEDURE ezcdf_close
   END INTERFACE

CONTAINS

   subroutine ezcdf_open(ncid, filename, opt, ier)
      ! Create/open cdf file
      ! 03/09/99 C.Ludescher
      !
      !  opt values:
      !    "R" -- readonly, existing file
      !    "W" -- write new file
      !    "M" -- "Modify" -- read/write data in existing file
      !                       (but structure of file does not change)
      !    "A" -- "Append" -- add new structure (define new data items)
      !                       in existing file.  Existing items can also
      !                       be read/written, but only after *new* items
      !                       are first defined and then written.
      !
      !  for both "W" and "A" modes, the file is opened in "define data mode".
      !
      include "netcdf.inc"
      INTEGER, intent(out) :: ncid
      character*(*), intent(in) :: filename
      character*1, intent(in) :: opt
      integer, optional, intent(out) :: ier
      integer :: status

      if (opt == 'w' .or. opt == 'W') then
         ! New file... start in "define mode".
         status = nf_create(filename, IOR(NF_CLOBBER, NF_64BIT_OFFSET), ncid)
         call handle_err(status, filename, 'cdfcrt', 'nf_create')
      else if (opt == 'm' .or. opt == 'M') then
         ! Open existing file for read/write modifications...
         status = nf_open(filename, nf_write, ncid)
         call handle_err(status, filename, 'cdfopn', 'nf_open')
      else if (opt == 'a' .or. opt == 'A') then
         ! Open existing file for read/write modifications...
         status = nf_open(filename, nf_write, ncid)
         call handle_err(status, filename, 'cdfopn', 'nf_open')
         if (status == NF_NOERR) then
            status = nf_redef(ncid)  ! start in "define mode".
            call handle_err(status, filename, 'cdfopn', 'nf_redef')
         end if
      else
         ! Open for readonly
         status = nf_open(filename, nf_nowrite, ncid)
         call handle_err(status, filename, 'cdfopn', 'nf_open')
      end if
      if (PRESENT(ier)) then
         if (status /= NF_NOERR) then
            ier = 1
         else
            ier = 0
         end if
      end if
      return
   end subroutine ezcdf_open

   subroutine ezcdf_close(ncid, ier)
      include "netcdf.inc"
      INTEGER, INTENT(in) ::  ncid
      integer, optional, intent(out) :: ier
      INTEGER status
      status = nf_close(ncid)
      call handle_err(status, ' ', 'cdfcls', 'nf_close')
      if (PRESENT(ier)) ier = status
   end subroutine ezcdf_close
!DEC$ ENDIF
END MODULE ezcdf_opncls
