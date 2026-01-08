! ================================================================================================================================================================================= !
! ---------------------------------------------------------------- Diagnostics for reading NEO's data on stella grids. ------------------------------------------------------------ !
! ================================================================================================================================================================================= ! 

module neoclassical_diagnostics
    implicit none

    ! Load debug flags?

    ! Make routines available to other modules.
    public :: write_neo_h_hat_on_stella_z_grid_diagnostic, write_neo_h_vmu_diagnostic, write_neo_phi_on_stella_z_grid_diagnostic
    public :: write_dneo_h_dz_diagnostic, write_dneo_phi_dz_diagnostic
    public :: write_wpol_diagnostic
    public :: write_dchidx_diagnostic_in_advance_wpol_routine

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
! ------------------------------------------------------------------------ Reads distributed dneo_h_dz data. ---------------------------------------------------------------------- !
! ================================================================================================================================================================================= !

    subroutine write_dneo_h_dz_diagnostic(dneo_h_dz, surface_index)
        use grids_z, only: nzgrid, zed
        use parallelisation_layouts, only: vmu_lo, iv_idx, imu_idx, is_idx
        use mp, only: proc0

        implicit none

        real, intent(in) :: dneo_h_dz(-nzgrid:nzgrid, vmu_lo%llim_proc:vmu_lo%ulim_proc, surface_index)
        integer, intent(in) :: surface_index

        integer :: unit
        character(len=256) :: filename
        integer :: iz, ivmu, iv, imu, is

        ! Only rank-0 writes diagnostics.
        if (.not. proc0) return
            unit = 99

            write(filename,'("dneo_h_dz_data",".dat")')
            

            open(unit=unit, file=filename, status='replace', action='write')

            write(unit,'(A)') '# =============================================================='
            write(unit,'(A)') '# Diagnostic output: dneo_h_dz (distributed vmu)'
            write(unit,'(A,I0)') '# surface_index = ', surface_index
            write(unit,'(A)') '# Columns:' 
            write(unit,'(A)') '#   iz   : z-grid index'
            write(unit,'(A)') '#   ivmu : local velocity-space index'
            write(unit,'(A)') '#   iv   : v_parallel index'
            write(unit,'(A)') '#   imu  : mu index'
            write(unit,'(A)') '#   is   : species index'
            write(unit,'(A)') '#   z    : stella z coordinate'
            write(unit,'(A)') '#   dneo_h_dz: neoclassical correction z derivative'
            write(unit,'(A)') '#'
            write(unit,'(A)') '# iz   ivmu   iv   imu   is        z              dneo_h_dz'
            write(unit,'(A)') '# ============================================================='

            do ivmu = vmu_lo%llim_proc, vmu_lo%ulim_proc
                iv  = iv_idx (vmu_lo, ivmu)
                imu = imu_idx(vmu_lo, ivmu)
                is  = is_idx (vmu_lo, ivmu)

                do iz = -nzgrid, nzgrid
                    write(unit,'(I6,1X,I6,1X,I4,1X,I4,1X,I4,1X,ES16.8,1X,ES16.8)') iz, ivmu, iv, imu, is, zed(iz), dneo_h_dz(iz, ivmu, 1)
                end do
            end do
            close(unit)
    end subroutine write_dneo_h_dz_diagnostic


! ================================================================================================================================================================================= !
! -------------------------------------------------------------------------------- Reads dneo_phi_dz data. ------------------------------------------------------------------------ !
! ================================================================================================================================================================================= !

    subroutine write_dneo_phi_dz_diagnostic(dneo_phi_dz, surface_index)
        use grids_z, only: nzgrid, zed
        use mp, only: proc0

        implicit none

        real, intent(in) :: dneo_phi_dz(-nzgrid:nzgrid, surface_index)
        integer, intent(in) :: surface_index

        integer :: unit
        character(len=256) :: filename
        integer :: iz

        ! Only rank-0 writes diagnostics.
        if (.not. proc0) return
            unit = 99

            write(filename,'("dneo_phi_dz_data",".dat")')
            
            open(unit=unit, file=filename, status='replace', action='write')

            write(unit,'(A)') '# ======================================================='
            write(unit,'(A)') '# Diagnostic output: dneo_phi_dz'
            write(unit,'(A,I0)') '# surface_index = ', surface_index
            write(unit,'(A)') '# Columns:' 
            write(unit,'(A)') '#   iz   : z-grid index'
            write(unit,'(A)') '#   z    : stella z coordinate'
            write(unit,'(A)') '#   dneo_phi_dz: neoclassical correction z derivative'
            write(unit,'(A)') '#'
            write(unit,'(A)') '# iz    z         dneo_phi_dz'
            write(unit,'(A)') '# ======================================================='

            do iz = -nzgrid, nzgrid
                write(unit,'(I6,1X,ES16.8,1X,ES16.8)') iz, zed(iz), dneo_phi_dz(iz, surface_index)
            end do
          
            close(unit)
    end subroutine write_dneo_phi_dz_diagnostic


! ================================================================================================================================================================================= !
! ---------------------------------------------------------------------- Reads wpol data in the init_wpol routine. ---------------------------------------------------------------- !
! ================================================================================================================================================================================= !

    subroutine write_wpol_diagnostic(wpol)
        use grids_z, only: nzgrid
        use grids_kxky, only: nalpha
        use parallelisation_layouts, only: vmu_lo, iv_idx, imu_idx, is_idx
        use mp, only: proc0

        implicit none

        real, intent(in) :: wpol(nalpha, -nzgrid:nzgrid, vmu_lo%llim_proc:vmu_lo%ulim_alloc)

        integer :: unit
        character(len=256) :: filename
        integer :: iz, ivmu, iv, imu, is, ia

        ! Only rank 0 writes out. 
        if (.not. proc0) return

        unit = 99
        filename = 'wpol_data.dat'

        open(unit=unit, file=filename, status='replace', action='write')

        ! Header
        write(unit,'(A)') '# ================================================='
        write(unit,'(A)') '# Diagnostic output: wpol'
        write(unit,'(A)') '# Columns:' 
        write(unit,'(A)') '#   iz    : z-grid index'
        write(unit,'(A)') '#   ivmu  : local vmu index'
        write(unit,'(A)') '#   iv    : v_parallel index'
        write(unit,'(A)') '#   imu   : mu index'
        write(unit,'(A)') '#   is    : species index'
        write(unit,'(A)') '#   wpol  : value'
        write(unit,'(A)') '# ================================================='

        ! Loop over all indices explicitly.
        
        do ivmu = vmu_lo%llim_proc, vmu_lo%ulim_proc
            iv  = iv_idx(vmu_lo, ivmu)
            imu = imu_idx(vmu_lo, ivmu)
            is  = is_idx(vmu_lo, ivmu)

            do iz = -nzgrid, nzgrid
                do ia = 1, nalpha
                    write(unit,'(I6,1X,I6,1X,I6,1X,I6,1X,I6,1X,ES16.8)') iz, ivmu, iv, imu, is, wpol(ia, iz, ivmu)
                end do
            end do
        end do

    close(unit)
    end subroutine write_wpol_diagnostic


! ================================================================================================================================================================================= !
! ----------------------------------------------------------- Reads dchidx data in the advance_wpol routine at a given timestep. -------------------------------------------------- !
! ================================================================================================================================================================================= !

    subroutine write_dchidx_diagnostic_in_advance_wpol_routine(dchidx)
        use grids_z, only: nzgrid, ntubes 
        use parallelisation_layouts, only: vmu_lo, iv_idx, imu_idx, is_idx
        use grids_kxky, only: nakx, naky
        use mp, only: proc0

        implicit none

        complex, intent(in) :: dchidx(:, :, -nzgrid:, :, vmu_lo%llim_proc:)

        integer :: unit
        character(len=256) :: filename
        integer :: iky, ikx, iz, it, ivmu
        integer :: iv, imu, is

        ! Only rank 0 writes out. 
        if (.not. proc0) return

        unit = 99        
        write(filename,'("dchidx_t.dat")')

        open(unit=unit, file=filename, status='replace', action='write')

        ! Header.
        write(unit,'(A)') '# ================================================='
        write(unit,'(A)') '# Diagnostic output: dchidx'
        write(unit,'(A)') '# Columns:' 
        write(unit,'(A)') '#   iy         : y mode number'
        write(unit,'(A)') '#   ix         : x mode number'
        write(unit,'(A)') '#   iz         : Local z index'
        write(unit,'(A)') '#   it         : flux tube index'
        write(unit,'(A)') '#   ivmu       : Local MPI index'
        write(unit,'(A)') '#   dchicdx_re : real value'
        write(unit,'(A)') '#   dchicdx_im : imaginary value'
        write(unit,'(A)') '# ================================================='

        ! Loop over all indices explicitly.
        
        do ivmu = vmu_lo%llim_proc, vmu_lo%ulim_proc
            iv  = iv_idx(vmu_lo, ivmu)
            imu = imu_idx(vmu_lo, ivmu)
            is  = is_idx(vmu_lo, ivmu)
  
            do iky = 1, naky
                do ikx = 1, nakx
                    do iz = -nzgrid, nzgrid
                        do it = 1, ntubes
                            write(unit,'(I6,1X,I6,1X,I6,1X,I6,1X,I6,1X,ES16.8,1X,ES16.8)') iky, ikx, iz, it, ivmu, real(dchidx(iky, ikx, iz, it, ivmu)), aimag(dchidx(iky, ikx, iz, it, ivmu))
                        end do
                    end do
                end do
            end do
        end do

    close(unit)
    end subroutine write_dchidx_diagnostic_in_advance_wpol_routine


! ================================================================================================================================================================================= !
! ----------------------------------------------------------------------------------- End Module. --------------------------------------------------------------------------------- !
! ================================================================================================================================================================================= !

end module neoclassical_diagnostics
