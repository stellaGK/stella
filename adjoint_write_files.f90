module adjoint_write_files

   implicit none
   public :: write_files_derivative
   public :: write_files_omega

   private
contains

   subroutine write_g_start
      use dist_fn_arrays, only: gnew

      use vpamu_grids, only: integrate_species
      use mp, only: sum_allreduce
      use volume_averages, only: fieldline_average

      use species, only: nspec
      use kt_grids, only: naky, nakx
      use zgrid, only: nzgrid, ntubes

      implicit none

      complex, dimension(:, :, :, :), allocatable :: int_tubes
      complex, dimension(:, :), allocatable :: gsum
      real, dimension(:), allocatable :: wgts
      integer :: iz, it

      allocate (int_tubes(naky, nakx, -nzgrid:nzgrid, ntubes)); int_tubes = 0.0
      allocate (gsum(naky, nakx)); gsum = 0.0
      allocate (wgts(nspec)); wgts = 1.0

      do it = 1, ntubes
         do iz = -nzgrid, nzgrid
            call integrate_species(gnew(:, :, iz, it, :), iz, wgts, int_tubes(:, :, iz, it), it, reduce_in=.false.)
         end do
      end do

      call sum_allreduce(int_tubes)

      call fieldline_average(int_tubes, gsum)

      open (14, file="adjoint_files/adjoint_ginit.dat", status="replace", action="write")
      write (14, *) abs(gsum(1, 1))
      close (14)

      deallocate (gsum, wgts, int_tubes)

   end subroutine write_g_start

   subroutine write_g_final

      use dist_fn_arrays, only: gnew

      use vpamu_grids, only: integrate_species
      use mp, only: sum_allreduce
      use volume_averages, only: fieldline_average

      use species, only: nspec
      use kt_grids, only: naky, nakx
      use zgrid, only: nzgrid, ntubes

      implicit none

      complex, dimension(:, :, :, :), allocatable :: int_tubes
      complex, dimension(:, :), allocatable :: gsum
      real, dimension(:), allocatable :: wgts
      integer :: iz, it

      allocate (int_tubes(naky, nakx, -nzgrid:nzgrid, ntubes)); int_tubes = 0.0
      allocate (gsum(naky, nakx)); gsum = 0.0
      allocate (wgts(nspec)); wgts = 1.0

      do it = 1, ntubes
         do iz = -nzgrid, nzgrid
            call integrate_species(gnew(:, :, iz, it, :), iz, wgts, int_tubes(:, :, iz, it), it, reduce_in=.false.)
         end do
      end do
      call sum_allreduce(int_tubes)

      call fieldline_average(int_tubes, gsum)

      open (15, file="adjoint_files/adjoint_gend.dat", status="replace", action="write")
      write (15, *) abs(gsum(1, 1))
      close (15)

      deallocate (gsum, wgts, int_tubes)

   end subroutine write_g_final

   subroutine write_final_time(istep_final)

!    use stella_time, only: code_dt
      use run_parameters, only: nstep
      use mp, only: proc0

      implicit none

      real, intent(in) :: istep_final

      if(proc0) then 
         open (16, file="adjoint_files/adjoint_final_time.dat", status="replace", action="write")
         write (16, *) istep_final - 1
         write (16, *) nstep
         close (16)
      end if

   end subroutine write_final_time

   subroutine write_files_omega(istep_final)

      use adjoint_field_arrays, only: omega_g
      use mp, only: proc0

      implicit none

      real, intent(in) :: istep_final

      if(proc0) then 
         open (13, file="adjoint_files/adjoint_omega.dat", status="replace", action="write")
         write (13, *) real(omega_g)
         close (13)
      end if

      call write_g_start
      call write_final_time(istep_final)

   end subroutine write_files_omega

   subroutine write_files_derivative(adjoint_var, derivative, new_file)

      use adjoint_field_arrays, only: omega_g
      use stella_geometry, only: geo_surf
      use mp, only: proc0

      implicit none

      integer, intent(in) :: adjoint_var
      complex, dimension(:, :), intent(in) :: derivative
      logical, intent(in) :: new_file

      if (new_file) then
         if(proc0) then 
            open (12, file="adjoint_files/adjoint_derivatives.dat", status="replace", action="write")!,position="replace")
            write (12, *) real(derivative(1, 1))
            close (12)
         end if
      else
         if(proc0) then 
            open (12, file="adjoint_files/adjoint_derivatives.dat", status="unknown", action="write", position="append")
            write (12, *) real(derivative(1, 1))
            close (12)
         end if
      end if

      if (.not. new_file) call write_g_final

   end subroutine write_files_derivative

end module adjoint_write_files
