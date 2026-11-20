!###############################################################################
!                                    TIMERS                                    
!###############################################################################
! This module gives Fortran access to the named timers of the C library <region_timer>
! (STELLA_CODE/region_timer). Regions are created on their first <region_start> call.
!###############################################################################
module timers

   use iso_c_binding, only: c_char, c_null_char, c_double

   implicit none
   
   ! Named timers from the C library <region_timer> (STELLA_CODE/region_timer)
   public :: region_start
   public :: region_end
   public :: region_total_time
   
   private

   !----------------------------------------------------------------------------
   !                  Interface to the C library <region_timer>                  
   !----------------------------------------------------------------------------
   ! Usage:
   !    call region_start('fields')
   !    ...
   !    call region_end('fields')
   !    t = region_total_time('fields')   ! accumulated time in seconds
   ! The strings are converted to null-terminated C strings in the wrappers below.
   
   interface
      subroutine c_region_start(name) bind(C, name="region_start")
         import :: c_char
         character(kind=c_char), dimension(*), intent(in) :: name
      end subroutine c_region_start
      
      subroutine c_region_end(name) bind(C, name="region_end")
         import :: c_char
         character(kind=c_char), dimension(*), intent(in) :: name
      end subroutine c_region_end
      
      function c_region_total_time(name) bind(C, name="region_total_time") result(total_time)
         import :: c_char, c_double
         character(kind=c_char), dimension(*), intent(in) :: name
         real(c_double) :: total_time
      end function c_region_total_time
   end interface

contains

   subroutine region_start(name)
      character(*), intent(in) :: name
      call c_region_start(trim(name)//c_null_char)
   end subroutine region_start

   subroutine region_end(name)
      character(*), intent(in) :: name
      call c_region_end(trim(name)//c_null_char)
   end subroutine region_end

   ! Total time (in seconds) spent in a region, 0 if the region was never entered
   function region_total_time(name) result(total_time)
      character(*), intent(in) :: name
      real :: total_time
      total_time = real(c_region_total_time(trim(name)//c_null_char))
   end function region_total_time

end module timers
