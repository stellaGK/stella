! ================================================================================================================================================================================= !
! ---------------------------------------------------------------- Diagnostics for reading NEO's data on stella grids. ------------------------------------------------------------ !
! ================================================================================================================================================================================= ! 

module neoclassical_diagnostics
    implicit none

    ! Load debug flags?

    ! Make routines available to other modules.
    public :: write_neo_h_hat_on_stella_z_grid_diagnostic, write_neo_h_vmu_diagnostic, write_neo_phi_on_stella_z_grid_diagnostic

contains 

! ================================================================================================================================================================================= !
! ---------------------------------------------- Reads NEO h hat data on stellas z grid. The data is still on the NEO velocity grids. --------------------------------------------- !
! ================================================================================================================================================================================= !

    subroutine write_neo_h_hat_on_stella_z_grid_diagnostic(neo_h_hat_in, neo_grid, surface_index, neo_h_hat_z_grid, suffix)
        use grids_z, only: nzgrid, zed
        use NEO_interface, only: neo_grid_data
    
        implicit none

        real, intent(in)  :: neo_h_hat_in(:, :, :, :, :)    
        real, intent(out) :: neo_h_hat_z_grid(-nzgrid:, :, :, :, :)

        type(neo_grid_data), intent(in) :: neo_grid
        integer, intent(in) :: surface_index
        integer :: ix, ie, is, iz                                 
        character(len=*), intent(in), optional :: suffix              
        integer :: unit
        character(len=256) :: filename

        unit = 99 

        if (present(suffix)) then
            write(filename,'("neo_h_hat_z_grid_",A,"_surf_",I0,".dat")') trim(suffix), surface_index
        else
            write(filename,'("neo_h_hat_z_grid_surf_",I0,".dat")') surface_index
        end if

        open(unit=unit, file=filename, status='replace', action='write')

        write(unit,'(A)') '# Diagnostic output: neo_h_hat_z_grid'
        write(unit,'(A,I0)') '# surface_index = ', surface_index
        write(unit,'(A)') '# Columns:'
        write(unit,'(A)') '#   iz   : z-grid index'
        write(unit,'(A)') '#   ix   : NEO xi index'
        write(unit,'(A)') '#   ie   : NEO energy index'
        write(unit,'(A)') '#   is   : species index'
        write(unit,'(A)') '#   z    : stella z coordinate'
        write(unit,'(A)') '#   hhat : interpolated neo_h_hat on stella z-grid'
        write(unit,'(A)') '#'
        write(unit,'(A)') '# iz   il   ie   is        z              hhat'
        ! -----------------------------------

        do ix = 1, neo_grid%n_xi + 1
            do ie = 1, neo_grid%n_energy + 1
                do is = 1, neo_grid%n_species
                    do iz = -nzgrid, nzgrid
                        write(unit,'(I6,1X,I4,1X,I4,1X,I4,1X,ES16.8,1X,ES16.8)') iz, ix, ie, is, zed(iz), neo_h_hat_z_grid(iz, ix, ie, is, surface_index)
                    end do
                end do
            end do
        end do
        close(unit)
    end subroutine write_neo_h_hat_on_stella_z_grid_diagnostic


! ================================================================================================================================================================================= !
! --------------------------------------------------------------- Reads NEO h hat data on stellas z and velocity grids. ----------------------------------------------------------- !
! ================================================================================================================================================================================= !

    subroutine write_neo_h_vmu_diagnostic(neo_h, surface_index, suffix)
        use grids_z, only: nzgrid, zed
        use parallelisation_layouts, only: vmu_lo, iv_idx, imu_idx, is_idx
        use mp, only: proc0

        implicit none

        real, intent(in) :: neo_h(-nzgrid:nzgrid, vmu_lo%llim_proc:vmu_lo%ulim_proc, surface_index)
        integer, intent(in) :: surface_index
        character(len=*), intent(in), optional :: suffix

        integer :: unit
        character(len=256) :: filename
        integer :: iz, ivmu, iv, imu, is

        ! Only rank-0 writes diagnostics.
        if (.not. proc0) return
            unit = 99

            if (present(suffix)) then
                write(filename,'("neo_h_vmu_",A,"_surf_",I0,".dat")') trim(suffix), surface_index
            else
                write(filename,'("neo_h_vmu_surf_",I0,".dat")') surface_index
            end if

            open(unit=unit, file=filename, status='replace', action='write')

            write(unit,'(A)') '# ================================================='
            write(unit,'(A)') '# Diagnostic output: neo_h (distributed vmu)'
            write(unit,'(A,I0)') '# surface_index = ', surface_index
            write(unit,'(A)') '# Columns:' 
            write(unit,'(A)') '#   iz   : z-grid index'
            write(unit,'(A)') '#   ivmu : local velocity-space index'
            write(unit,'(A)') '#   iv   : v_parallel index'
            write(unit,'(A)') '#   imu  : mu index'
            write(unit,'(A)') '#   is   : species index'
            write(unit,'(A)') '#   z    : stella z coordinate'
            write(unit,'(A)') '#   neo_h: neoclassical correction'
            write(unit,'(A)') '#'
            write(unit,'(A)') '# iz   ivmu   iv   imu   is        z              neo_h'
            write(unit,'(A)') '# -------------------------------------------------'

            do ivmu = vmu_lo%llim_proc, vmu_lo%ulim_proc
                iv  = iv_idx (vmu_lo, ivmu)
                imu = imu_idx(vmu_lo, ivmu)
                is  = is_idx (vmu_lo, ivmu)

                do iz = -nzgrid, nzgrid
                    write(unit,'(I6,1X,I6,1X,I4,1X,I4,1X,I4,1X,ES16.8,1X,ES16.8)') iz, ivmu, iv, imu, is, zed(iz), neo_h(iz, ivmu, 1)
                end do
            end do
            close(unit)
    end subroutine write_neo_h_vmu_diagnostic


! ================================================================================================================================================================================= !
! --------------------------------------------------------------- Reads NEO ϕ^1_0 data on stellas z and velocity grids. ----------------------------------------------------------- !
! ================================================================================================================================================================================= !

    subroutine write_neo_phi_on_stella_z_grid_diagnostic(neo_phi_in, neo_grid, surface_index, neo_phi)
        use grids_z, only: nzgrid, zed
        use NEO_interface, only: neo_grid_data
      
        implicit none

        real, intent(in)  :: neo_phi_in(:, :)    
        type(neo_grid_data), intent(in) :: neo_grid
        integer, intent(in) :: surface_index
        real, intent(out) :: neo_phi(-nzgrid:, :)        
        integer :: iz
         
        integer :: unit
        character(len=256) :: filename

        unit = 99 

        write(filename,'("neo_phi",I0,".dat")') surface_index
      
        open(unit=unit, file=filename, status='replace', action='write')

        write(unit,'(A)') '# Diagnostic output: neo_phi'
        write(unit,'(A,I0)') '# surface_index = ', surface_index
        write(unit,'(A)') '# Columns:'
        write(unit,'(A)') '#   iz   : z-grid index'
        write(unit,'(A)') '#   z    : stella z coordinate'
        write(unit,'(A)') '#   phi^1_0 : interpolated neo_phi on stella z-grid'
        write(unit,'(A)') '#'
        write(unit,'(A)') '# iz            z              hhat'
        ! -----------------------------------

        do iz = -nzgrid, nzgrid
            write(unit,'(I6,ES16.8,1X,ES16.8)') iz, zed(iz), neo_phi(iz, surface_index)
        end do
        close(unit)
    end subroutine write_neo_phi_on_stella_z_grid_diagnostic


! ================================================================================================================================================================================= !
! ----------------------------------------------------------------------------------- End Module. --------------------------------------------------------------------------------- !
! ================================================================================================================================================================================= !


end module neoclassical_diagnostics
