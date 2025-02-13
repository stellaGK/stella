!###############################################################################
!################################ MAIN PROGRAM #################################
!###############################################################################

! We define the main code as a module, since we want to use nested 
! subroutines, which is not allowed in a program.
program update_input_file_program
   use input_file, only: update_input_file
   implicit none   
   call update_input_file() 
end program update_input_file_program
