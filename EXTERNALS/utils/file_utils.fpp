# include "define.inc"

module file_utils

   implicit none

   private

   public :: init_file_utils
   ! subroutine init_file_utils (list, input, error, trin_run, name)
   ! logical, intent (out) :: list
   ! logical, intent (in), optional :: input, error, trin_run
   ! character(*), intent (in), optional :: name
   !   default: INPUT=.true., ERROR=.true., TRIN_RUN=.false., NAME="unknown"
   !   Set up run_name(s) and list_name for output files
   !   Open input file and strip comments, unless disabled with INPUT=.false.
   !   Open error output file, unless disabled with ERROR=.false.

   public :: init_job_name
   ! subroutine ...

   public :: finish_file_utils
   ! subroutine finish_file_utils
   !   Clean up files opened in init

   public :: run_name
   ! character(500) :: run_name
   !    Label for the run, taken from the command line

   public :: list_name
   ! character(500) :: list_name
   !    Label for the list, taken from the command line

   public :: input_unit
   ! function input_unit (nml)
   ! character(*), intent (in) :: nml
   ! integer :: input_unit
   !    Rewind the input file to start of namelist NML,
   !    and return its unit number

   public :: input_unit_exist
   ! function input_unit_exist (nml,exist)
   ! character(*), intent (in) :: nml
   ! integer :: input_unit
   !    Rewind the input file to start of namelist NML,
   !    and return its unit number, setexist=.true.
   !    If the namelist NML is not found, set exist=.false.

   public :: open_error_file
   public :: write_clean_input_file

   public :: error_unit
   ! function error_unit ()
   ! integer :: error_unit
   !    Return the error unit number

   public :: get_input_unit

   public :: open_output_file
   ! subroutine open_output_file (unit, ext)
   ! integer, intent (out) :: unit
   ! character (*), intent (in) :: ext
   !    Open a file with name made from the run_name with the EXT appended
   !    and return its unit number in UNIT

   public :: close_output_file
   ! subroutine close_output_file (unit)
   ! integer, intent (in) :: unit
   !    Close the file associated with UNIT from open_output_file

   public :: flush_output_file
   ! subroutine flush_output_file (unit)
   ! integer, intent (in) :: unit
   !    Close/open-append the file associated with UNIT from open_output_file

   public :: get_unused_unit
   ! subroutine get_unused_unit (unit)
   ! integer, intent (out) :: unit
   !    Return a unit number not associated with any file

   public :: get_indexed_namelist_unit
   ! subroutine get_indexed_namelist_unit (unit, nml, index)
   ! integer, intent (out) :: unit
   ! character (*), intent (in) :: nml
   ! integer, intent (in) :: index
   !    Copy namelist, NML // '_' // INDEX, from the input file to
   !    namelist, NML, in a temporary file, UNIT
   
   ! Allow us to switch input files
   public :: open_other_input_file

!  public :: num_input_lines

   public :: stdout_unit

   public :: runtype_option_switch
   public :: runtype_standalone
   public :: runtype_trinity
   public :: runtype_list
   public :: runtype_multibox

   character(500), pointer :: run_name
   character(500), target :: arun_name, job_name
   character(500) :: list_name
   integer, parameter :: stdout_unit = 6
   integer :: runtype_option_switch
   integer, parameter :: runtype_standalone = 0, &
                         runtype_list = 1, &
                         runtype_trinity = 2, &
                         runtype_multibox = 3

   integer, save :: input_unit_no, error_unit_no = stdout_unit
   integer, save, public :: num_input_lines

contains

   !============================================================================
   !======================== Initialize the file utils =========================
   !============================================================================  
   ! Read the input file name from the command line, and save it as <arun_name>.
   ! Save the name of the input file without it's extension ".in" as <run_name>.
   ! If the extension of the input file is ".list" or ".multi" set list = .true.
   
   ! Find out the [[run_name]], and use the run name to determine whether
   ! this is a [[list]] run (i.e. a list of runs has been given) or a [[Trinity]] run.
   ! If not, open the error file and call write_clean_input_file
   subroutine init_file_utils(list)
   
      implicit none
      
      logical, intent(out) :: list

      ! Get the name of the input file (arun_name) from the command line
      ! and set list = .true. if the input file name ends in ".list" or ".multi"
      call get_name_input_file(list)

      ! If the input file name ends in ".list" or ".multi"
      if (list) then
         list_name = arun_name

      ! If the input file name ends in ".in"
      else
      
         ! Save the name of the input file without it's extension ".in" as <run_name>. 
         call get_run_name()
         
         ! Open the error file <run_name>.error and set its unit number to <error_unit_no>. 
         call open_error_file()
         
         ! Read the user specified input file, and clean it up by removing
         ! comments from the file and by reading in nested input files.
         ! stella will only use the cleaned up file, not the file from the user. 
         call write_clean_input_file()
         
      end if

   end subroutine init_file_utils
   
   
!###############################################################################
!######################### NAME INPUT FILE (RUN_NAME) ##########################
!###############################################################################

   !============================================================================
   !==================== Name of the input file (arun_name) ====================
   !============================================================================ 
   subroutine get_name_input_file(list)
      ! This determines the type of run, by reading the name of the input file
      ! on the command line into [[arun_name]], and then looking at the extension. If
      ! the extension is .list, then [[list]] is set to .true.).

      use command_line, only: cl_getarg, cl_iargc 

      implicit none
      logical, intent(out) :: list
      integer :: l, ierr

      ! Initialize
      list = .false.
      
      ! Get the first argument from the command line and put it in <arun_name>
      if (cl_iargc() /= 0) then
         call cl_getarg(1, arun_name, l, ierr)
         if (ierr /= 0) then
            print *, "Error getting run name."
         end if
      end if
      
      ! Exit the program if no input file has been specified.
      if (l < 2 .or. l > 500) then 
         write(*,*) ' '; write(*,*) 'ERROR: Please specify an input file. For example:'; 
         write(*,*) '       >> mpirun -np 2 stella input.in '; write(*,*) ' '; stop
      end if

      ! Check if <arun_name> end in ".list"
      if (l > 5 .and. arun_name(l - 4:l) == ".list") then
         list = .true.
         runtype_option_switch = runtype_list
      end if

      ! Check if <arun_name> end in ".multi"
      if (l > 6 .and. arun_name(l - 5:l) == ".multi") then
         list = .true.
         runtype_option_switch = runtype_multibox
      end if

   end subroutine get_name_input_file

   ! Remove the extension from <arun_name> and put it in <run_name>
   subroutine get_run_name
   
      implicit none
      
      integer :: l

      l = len_trim(arun_name)
      if (l > 3 .and. arun_name(l - 2:l) == ".in") then
         arun_name = arun_name(1:l - 3)
      end if
      run_name => arun_name

   end subroutine get_run_name

   ! Define the <run_name> or <job_name> on all processors.
   ! If a list of input files is used, the jobs are spread out over the
   ! processors, and they can have different values of <job_name>.
   subroutine init_job_name(jobname)
   
      implicit none
      
      character(len=500), intent(in) :: jobname
      
      job_name = trim(jobname)
      run_name => job_name
      
   end subroutine init_job_name
   
!###############################################################################
!############################# CLEANED INPUT FILE ##############################
!###############################################################################

   !============================================================================
   !========================== Create cleaned input file =======================
   !============================================================================
   ! Open the input file <run_name>.in, strip out any comments, and write the  
   ! resulting lines into the file .<run_name>.in. 
   ! 
   ! Note that in the input file we can include the line:
   ! 
   ! !include other_input_file.in
   ! 
   ! which will allow us to split up the input file in smaller input files.
   !============================================================================
   subroutine write_clean_input_file()

      implicit none
      
      character(500) :: line
      
      ! To hold position of slash in input_file_name
      integer :: ind_slash
     
      ! In <stack> we will save the unit number of the additional input
      ! files which are mentioned through !include other_input_file.in
      ! Only allow up to 10 included/nested input files
      integer, parameter :: stack_size = 10
      integer, dimension(stack_size) :: stack
      integer :: stack_ptr
      
      ! The user specified input file has unit number <in_unit>
      ! The cleaned up input file has unit number <out_unit> = <input_unit_no>
      integer :: in_unit, out_unit, iostat
      
      !-------------------------------------------------------------------------

      ! Get a unit number for the "input_file_name.in" that the user wrote
      call get_unused_unit(in_unit)
      
      ! Open <run_name>.in
      open (unit=in_unit, file=trim(run_name)//".in", status="old", &
            action="read", iostat=iostat)
      if (iostat /= 0) then
         write(*,"(a)") "Could not open input file: "//trim(run_name)//".in"
      end if

      ! Get a unit number for the cleaned .<run_name>.in file we will write
      call get_unused_unit(out_unit)

      ! Determine if '/' is present in <run_name> and if so what position, 
      ! e.g. 'folder/<run_name>.in', will be split into path_to_file and file
      ind_slash = index(run_name, "/", .True.)
      
      ! No slash in <run_name>
      if (ind_slash == 0) then  
         open (unit=out_unit, file="."//trim(run_name)//".in")
         
      ! Slash in <run_name>
      else
         open (unit=out_unit, file=trim(run_name(1:ind_slash))//"."//trim(run_name(ind_slash + 1:))//".in")
      end if

      ! Initialize
      iostat = 0
      stack_ptr = 0
      num_input_lines = 0
      
      ! Read the input files (or nested input files) one line at a time
      ! The input file can be split into max 10 smaller input files
      ! which are included through '!include <input_file_name_small.in>'
      do
         read (unit=in_unit, fmt="(a)", iostat=iostat) line
         
         ! If we reached the end of the file, then iostat will not be 0,
         ! then change the user_input_unit_number to the previous one in the 
         ! <stack> list, which contains a list of input files that need to  
         ! be included. cycle will restart the do loop, exit will stop it
         if (iostat /= 0) then
            if (stack_ptr <= 0) exit
            close (unit=in_unit)
            iostat = 0
            in_unit = stack(stack_ptr)
            stack_ptr = stack_ptr - 1
            cycle
         end if
          
         ! We can split the input file up in upto 10 little input files.
         ! In the main input file we include the other files through:
         ! !include <input_file_name_small.in>
         if (line(1:9) == "!include ") then
         
            ! If more than 10 !include statements are found, stop reading them
            if (stack_ptr >= stack_size) then
               write(*,"(a)") "The !include statement at the start of the file is ignored"
               write(*,"(a)") "because the nesting is too deep (use max 10 files). " 
               write(*,"(a)") "Ignored file: "//trim(line)
               cycle
            end if
            
            ! Save the unit number of the current input file in <stack>, and 
            ! start reading the file mentioned in !include <input_file_name_small.in>
            ! To get the name of the file, remove '!include' with trim(line(10:))
            stack_ptr = stack_ptr + 1
            stack(stack_ptr) = in_unit
            call get_unused_unit(in_unit)
            open (unit=in_unit, file=trim(line(10:)), status="old", &
                  action="read", iostat=iostat)
            if (iostat /= 0) then
               write(*,"(a)") "The !include statement at the start of the file is ignored"
               write(*,"(a)") "because the following input file could not be read: " 
               write(*,"(a)") "     "//trim(line)
               in_unit = stack(stack_ptr)
               stack_ptr = stack_ptr - 1
               cycle
            end if
            cycle
         end if
         
         ! Remove comments from the file
         call strip_comments(line)
         
         ! Write this line to the .<run_name>.in file 
         write (unit=out_unit, fmt="(a)") trim(line)
         num_input_lines = num_input_lines + 1
         
      end do
      
      ! Close the original input file, we will not use it anymore
      ! Since stella will only read the cleaned up input file 
      close (unit=in_unit)

      ! Save the unit number of the cleaned input file to <input_unit_no>
      input_unit_no = out_unit
      
   end subroutine write_clean_input_file
   
   ! Notice that in <write_clean_input_file> we never closed the input file
   ! since we will read it throughout stella, and we save its unit number
   ! under <input_unit_no>. We now want to switch which input file is being read.
   subroutine open_other_input_file(path_input_file)
   
      implicit none 
      
      character(500), intent(in) :: path_input_file
      integer :: iostat 
      
      ! Close the current input file which is open
      close (unit=input_unit_no) 
      
      ! Get a unit number for <path_input_file>
      call get_unused_unit(input_unit_no)
      
      ! Open the new input file
      open (unit=input_unit_no, file=path_input_file, status="old", action="read", iostat=iostat)
      if (iostat /= 0) then
         write(*,*) "Could not open switched input file: "
         write(*,*) "    ", path_input_file
      end if
   
   end subroutine open_other_input_file
   
!###############################################################################
!############################## OPEN/CLOSE FILES ###############################
!###############################################################################

   ! Get an unused unit number for I/O.
   subroutine get_unused_unit(unit) 
      implicit none
      integer, intent(out) :: unit
      logical :: od
      unit = 50
      do
         inquire (unit=unit, opened=od)
         if (.not. od) return
         unit = unit + 1
      end do
   end subroutine get_unused_unit

   !==============================================
   !============= OPEN OUTPUT FILE ===============
   !==============================================
   ! Open an output file to write data (replacing or appending any existing)
   ! The name is [[run_name]] + [[ext]], and set [[unit]] to the
   ! unit number of that output file.
   subroutine open_output_file(unit, ext, overwrite_in)

      implicit none

      integer, intent(out) :: unit
      logical, intent(in), optional :: overwrite_in
      logical :: overwrite
      character(*), intent(in) :: ext
      character(500) :: hack

      ! Initiate the optional argument
      if (present(overwrite_in)) then
         overwrite = overwrite_in
      else
         overwrite = .true.
      end if

      ! Get a unit for the output file that is not currently in use
      call get_unused_unit(unit)

      ! Create the name of the output file
      hack = trim(run_name)//ext

      ! If overwrite==True: Create a new output file or replace the existing file
      ! If overwrite==False: Append data to the already existing output file
      if (overwrite) then
         open (unit=unit, file=trim(hack), status="replace", action="write")
      else
         open (unit=unit, file=trim(hack), status="unknown", action="write", position="append")
      end if

   end subroutine open_output_file

   !==============================================
   !============= CLOSE OUTPUT FILE ==============
   !==============================================
   ! Close the output file identified by [[unit]].
   subroutine close_output_file(unit)
      implicit none
      integer, intent(in) :: unit
      close (unit=unit)
   end subroutine close_output_file

   subroutine flush_output_file(unit)
      implicit none
      integer, intent(in) :: unit
      character(len=500) :: fname
      inquire (unit, name=fname)
# if FCOMPILER == _XL_
      call flush_(unit)
# elif FCOMPILER == _NAG_
      close (unit=unit)
      open (unit=unit, file=trim(fname), status="old", action="write", position="append")
# else
      call flush (unit)
# endif
   end subroutine flush_output_file

   subroutine open_error_file()
      implicit none
      error_unit_no = 0
      if (run_name /= "unknown") then
         call open_output_file(error_unit_no, ".error")
      end if
   end subroutine open_error_file

   subroutine strip_comments(line)
      implicit none
      character(*), intent(in out) :: line
      logical :: in_single_quotes, in_double_quotes
      integer :: i, length

      length = len_trim(line)
      i = 1
      in_single_quotes = .false.
      in_double_quotes = .false.
      loop: do
         if (in_single_quotes) then
            if (line(i:i) == "'") in_single_quotes = .false.
         else if (in_double_quotes) then
            if (line(i:i) == '"') in_double_quotes = .false.
         else
            select case (line(i:i))
            case ("'")
               in_single_quotes = .true.
            case ('"')
               in_double_quotes = .true.
            case ("!")
               i = i - 1
               exit loop
            end select
         end if
         if (i >= length) exit loop
         i = i + 1
      end do loop
      line = line(1:i)
   end subroutine strip_comments


   subroutine finish_file_utils
      implicit none
      if (input_unit_no > 0) then
         close (unit=input_unit_no)
         input_unit_no = -1
      end if
      if (error_unit_no > 0 .and. error_unit_no /= 6) then
         close (unit=error_unit_no)
         error_unit_no = -1
      end if
   end subroutine finish_file_utils
   
!###############################################################################
!################################## NAMELISTS ##################################
!###############################################################################

   function input_unit(nml)
      implicit none
      character(*), intent(in) :: nml
      character(len(nml)) :: nml_upper
      integer :: input_unit, iostat
      character(500) :: line
      intrinsic adjustl, trim
      input_unit = input_unit_no
      nml_upper = str_to_upper_case(nml)
      if (input_unit_no > 0) then
         rewind (unit=input_unit_no)
         do
            read (unit=input_unit_no, fmt="(a)", iostat=iostat) line
            if (iostat /= 0) then
               rewind (unit=input_unit_no)
               exit
            end if
            if (trim(adjustl(line)) == "&"//nml) then
               backspace (unit=input_unit_no)
               return
            end if
            if (trim(adjustl(line)) == "&"//nml_upper) then
               backspace (unit=input_unit_no)
               return
            end if
         end do
      end if
      write (unit=error_unit_no, fmt="('Could not find namelist: ',a)") nml
      write (unit=*, fmt="('Could not find namelist: ',a)") nml
   end function input_unit

   function input_unit_exist(nml, exist)
      implicit none
      character(*), intent(in) :: nml
      logical, intent(out) :: exist
      integer :: input_unit_exist, iostat
      character(len(nml)) :: nml_upper
      character(500) :: line
      intrinsic adjustl, trim
      
      input_unit_exist = input_unit_no
      nml_upper = str_to_upper_case(nml)
      exist = .true.
      if (input_unit_no > 0) then
         rewind (unit=input_unit_no)
         do
            read (unit=input_unit_no, fmt="(a)", iostat=iostat) line
            if (iostat /= 0) then
               rewind (unit=input_unit_no)
               exit
            end if
            if (trim(adjustl(line)) == "&"//nml) then
               backspace (unit=input_unit_no)
               return
            end if
            if (trim(adjustl(line)) == "&"//nml_upper) then
               backspace (unit=input_unit_no)
               return
            end if
         end do
      end if
      exist = .false.
   end function input_unit_exist
   
   function str_to_upper_case(str) result(str_upper)

      implicit none

      character(*), intent(in) :: str
      character(len(str)) :: str_upper
      integer :: ic, i
      character(26), parameter :: cap = 'ABCDEFGHIJKLMNOPQRSTUVWXYZ_'
      character(26), parameter :: low = 'abcdefghijklmnopqrstuvwxyz_'

      str_upper = str
      do i = 1, len_trim(str)
        ic = index(low, str(i:i))
        if (ic > 0) str_upper(i:i) = cap(ic:ic)
      end do

   end function str_to_upper_case

   function error_unit()
      implicit none
      integer :: error_unit
      error_unit = error_unit_no
   end function error_unit

   subroutine get_input_unit(unit)
      implicit none
      integer, intent(out) :: unit

      unit = input_unit_no

   end subroutine get_input_unit

   subroutine get_indexed_namelist_unit(unit, nml, index_in)
      implicit none
      integer, intent(out) :: unit
      character(*), intent(in) :: nml
      integer, intent(in) :: index_in
      character(500) :: line
      integer :: iunit, iostat, in_file
      integer :: ind_slash
      logical :: exist

      call get_unused_unit(unit)
!    open (unit=unit, status="scratch", action="readwrite")

      !Determine if '/' is in input name and if so what position
      !in the string is the last one (i.e. split run_name into path_to_file and file)
      ind_slash = index(run_name, "/", .True.)
      if (ind_slash == 0) then !No slash in name
         !Original behaviour
         open (unit=unit, file="."//trim(run_name)//".scratch")
      else
         !General behaviour
         open (unit=unit, file=trim(run_name(1:ind_slash))//"."//trim(run_name(ind_slash + 1:))//".scratch")
      end if

      write (line, *) index_in
      line = nml//"_"//trim(adjustl(line))
      in_file = input_unit_exist(trim(line), exist)
      if (exist) then
         iunit = input_unit(trim(line))
      else
         write (6, *) "get_indexed_namelist: following namelist not found ", trim(line)
         return
      end if

      read (unit=iunit, fmt="(a)") line
      write (unit=unit, fmt="('&',a)") nml

      do
         read (unit=iunit, fmt="(a)", iostat=iostat) line
         if (iostat /= 0 .or. trim(adjustl(line)) == "/") exit
         write (unit=unit, fmt="(a)") trim(line)
      end do
      write (unit=unit, fmt="('/')")
      rewind (unit=unit)
   end subroutine get_indexed_namelist_unit

end module file_utils
