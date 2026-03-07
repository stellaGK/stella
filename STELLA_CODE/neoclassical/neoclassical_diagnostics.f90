! ================================================================================================================================================================================= !
! ---------------------------------------------------------------- Diagnostics for reading NEO's data on stella grids. ------------------------------------------------------------ !
! ================================================================================================================================================================================= ! 

module neoclassical_diagnostics
    implicit none

    ! Load debug flags?

    ! Make routines available to other modules.
    public :: write_neo_distribution_on_stella_grids_diagnostic, write_neo_phi_on_stella_z_grid_diagnostic

contains 

! ================================================================================================================================================================================= !
! --------------------------------------------------------------- Reads NEO h hat data on stellas z and velocity grids. ----------------------------------------------------------- !
! ================================================================================================================================================================================= !

    subroutine write_neo_distribution_on_stella_grids_diagnostic(neo_grid, neo_h_global, surface_index, file_tag)
        ! Grids.
        use grids_z, only: nzgrid, zed
        use grids_velocity, only: nvpa, nmu

        ! Parallelisation.
        use parallelisation_layouts, only: vmu_lo, iv_idx, imu_idx, is_idx

        ! MP. 
        use mp, only: proc0

        ! NEO.
        use NEO_interface, only: neo_grid_data

        implicit none

        type(neo_grid_data), intent(in) :: neo_grid
        real, intent(in) :: neo_h_global(-nzgrid:, :, :, :, :)
        integer, intent(in) :: surface_index
        character(len=*), intent(in) :: file_tag
        
        integer :: unit
        character(len=256) :: filename
        integer :: iz, iv, imu, is

        if (len_trim(file_tag) > 0) then
            filename = trim(file_tag) // ".dat"
        end if

        unit = 99

        open(newunit=unit, file=trim(filename), status='replace', action='write')

        write(unit,'(A)') '# ================================================='
        write(unit,'(A)') '# Diagnostic output: neo_h_global'
        write(unit,'(A,I0)') '# surface_index = ', surface_index
        write(unit,'(A)') '# Columns:' 
        write(unit,'(A)') '#   iz   : z-grid index'
        write(unit,'(A)') '#   iv   : v_parallel index'
        write(unit,'(A)') '#   imu  : mu index'
        write(unit,'(A)') '#   is   : species index'
        write(unit,'(A)') '#   z    : stella z coordinate'
        write(unit,'(A)') '#   neo_h_global: neoclassical correction'
        write(unit,'(A)') '#'
        write(unit,'(A)') '# iz   iv   imu   is        z                 neo_h'
        write(unit,'(A)') '# ================================================='

        do iz = -nzgrid, nzgrid
            do iv = 1, nvpa
                do imu = 1, nmu
                    do is = 1, neo_grid%n_species
                        write(unit,'(I6,1X,I4,1X,I4,1X,I4,1X,ES16.8,1X,ES16.8)') iz, iv, imu, is, zed(iz), neo_h_global(iz, iv, imu, is, 1)
                    end do
                end do
            end do
        end do
        
        close(unit)
    end subroutine write_neo_distribution_on_stella_grids_diagnostic


! ================================================================================================================================================================================= !
! --------------------------------------------------------------- Reads NEO ϕ^1_0 data on stellas z and velocity grids. ----------------------------------------------------------- !
! --------------------------------------------------------------- Will also take a ϕ^1_0 derivative as input if needed. ----------------------------------------------------------- !
! ================================================================================================================================================================================= !

    subroutine write_neo_phi_on_stella_z_grid_diagnostic(neo_grid, neo_phi, file_tag)
        use grids_z, only: nzgrid, zed
        use NEO_interface, only: neo_grid_data
      
        implicit none

        type(neo_grid_data), intent(in) :: neo_grid
        real, intent(in) :: neo_phi(-nzgrid:)        
        character(len=*), intent(in) :: file_tag 
        
        integer :: iz
        integer :: unit
        character(len=256) :: filename

        if (len_trim(file_tag) > 0) then
            filename = trim(file_tag) // ".dat"
        end if

        open(newunit=unit, file=trim(filename), status='replace', action='write')

        write(unit,'(A)') '# Diagnostic output: neo_phi'
        write(unit,'(A)') '# Columns:'
        write(unit,'(A)') '#    iz   : z-grid index'
        write(unit,'(A)') '#    z    : stella z coordinate'
        write(unit,'(A)') '#    phi^1_0 : interpolated neo_phi on stella z-grid'
        write(unit,'(A)') '#'
        write(unit,'(A)') '# iz            z                phi'
        
        do iz = -nzgrid, nzgrid
            write(unit,'(I6, 1X, ES16.8, 1X, ES16.8)') iz, zed(iz), neo_phi(iz)
        end do
        
        close(unit)
    end subroutine write_neo_phi_on_stella_z_grid_diagnostic


! ================================================================================================================================================================================= !
! -------------------------------------------------------- Reads the moments of the neoclassical H_1 distribution on stellas z grid. ---------------------------------------------- !
! ================================================================================================================================================================================= !

    subroutine write_distribution_moment_diagnostic(neo_grid, moment, file_tag)
        ! Grids.
        use grids_z, only: nzgrid, zed

        ! NEO. 
        use NEO_interface, only: neo_grid_data
      
        implicit none

        type(neo_grid_data), intent(in) :: neo_grid
        real, intent(in) :: moment(-nzgrid:, :)        
        character(len=*), intent(in) :: file_tag 
        
        integer :: iz, is, ia
        integer :: unit
        character(len=256) :: filename

        if (len_trim(file_tag) > 0) then
            filename = trim(file_tag) // ".dat"
        end if

        open(newunit=unit, file=trim(filename), status='replace', action='write')

        write(unit,'(A)') '# Diagnostic output: H_1 moment'
        write(unit,'(A)') '#  Columns:' 
        write(unit,'(A)') '#    is   : species index'
        write(unit,'(A)') '#    iz   : z-grid index'
        write(unit,'(A)') '#    z    : stella z coordinate'
        write(unit,'(A)') '#    mom  : H_1 velocity moment'
        write(unit,'(A)') '#'
        write(unit,'(A)') '#  is            iz            z                mom'
        
        do is = 1, neo_grid%n_species
            do iz = -nzgrid, nzgrid
                write(unit,'(I6, 1X, I6, 1X, ES16.8, 1X, ES16.8)') is, iz, zed(iz), moment(iz, is)
            end do
        end do
        
        close(unit)
    end subroutine write_distribution_moment_diagnostic


! ================================================================================================================================================================================= !
! ----------------------------------------------------------------------------------- End Module. --------------------------------------------------------------------------------- !
! ================================================================================================================================================================================= !

end module neoclassical_diagnostics
