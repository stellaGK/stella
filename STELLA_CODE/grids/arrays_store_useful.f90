module arrays_store_useful

    use mpi
    use stella_common_types, only: response_matrix_type

    implicit none
    
    ! Velocity-dependent coefficients used in the equations
    public :: kperp2, dkperp2dr
    public :: time_gke
    public :: time_parallel_nl

    public :: wdriftinit, wstarinit, parnlinit, &
            radialinit, driftimpinit

    public :: wstar, wstarp
    public :: wdriftx_g, wdrifty_g
    public :: wdriftx_phi, wdrifty_phi
    public :: wdriftx_bpar, wdrifty_bpar
    public :: wdriftpx_g, wdriftpy_g
    public :: wdriftpx_phi, wdriftpy_phi

    ! ! Arrays without velocity dependence. Used mostly in field calculations.
    public :: response_matrix, response_window
    public :: shift_state
    public :: gamtot, dgamtotdr
    public :: gamtot13, gamtot31, gamtot33
    public :: gamtotinv11, gamtotinv13, gamtotinv31, gamtotinv33
    public :: gamtot3
    public :: apar_denom
    public :: theta
    public :: c_mat
    public :: exclude_boundary_regions_qn
    public :: tcorr_source_qn, exp_fac_qn
    public :: qn_window, qn_zf_window
    public :: gamtot_h, gamtot3_h, efac, efacp
    public :: time_field_solve
    
    !----------------------------------------------------------------------------
    ! For GK eqn
    !----------------------------------------------------------------------------
    real, dimension(:, :, :), allocatable :: wstar, wstarp
    ! (nalpha, -nzgrid:nzgrid, -vmu-layout-)

    real, dimension(:, :, :), allocatable :: wdriftx_g, wdrifty_g
    real, dimension(:, :, :), allocatable :: wdriftx_phi, wdrifty_phi
    real, dimension(:, :, :), allocatable :: wdriftx_bpar, wdrifty_bpar

    real, dimension(:, :, :), allocatable :: wdriftpx_g, wdriftpy_g
    real, dimension(:, :, :), allocatable :: wdriftpx_phi, wdriftpy_phi
    ! (nalpha, -nzgrid:nzgrid, -vmu-layout-)

    !> dkperp2dr will contain the radial variation of kperp2
    real, dimension(:, :, :, :), allocatable :: kperp2, dkperp2dr

    !> for time advance
    real, dimension(2, 10) :: time_gke = 0.
    real, dimension(2, 2) :: time_parallel_nl = 0.

    logical :: wdriftinit, wstarinit, parnlinit, &
            radialinit, driftimpinit

    !----------------------------------------------------------------------------
    ! For field solves 
    !----------------------------------------------------------------------------
    type(response_matrix_type), dimension(:), allocatable :: response_matrix
    integer :: response_window = MPI_WIN_NULL
    real, dimension(:), allocatable :: shift_state
    real, dimension(:, :, :), allocatable :: gamtot, dgamtotdr
    real, dimension(:, :, :), allocatable :: gamtot13, gamtot31, gamtot33
    real, dimension(:, :, :), allocatable :: gamtotinv11, gamtotinv13, gamtotinv31, gamtotinv33
    real, dimension(:, :), allocatable :: gamtot3
    real, dimension(:, :, :), allocatable ::  apar_denom
    complex, dimension(:, :, :), allocatable :: theta
    ! (nakx, nakx, -nzgrid:nzgrid)
    complex, dimension(:, :), allocatable :: c_mat
    ! (nakx, nakx)
    !variables needed for the source
    logical :: exclude_boundary_regions_qn
    real :: tcorr_source_qn, exp_fac_qn
    integer :: qn_window = MPI_WIN_NULL, qn_zf_window = MPI_WIN_NULL
    real :: gamtot_h, gamtot3_h, efac, efacp
    real, dimension(2, 5) :: time_field_solve = 0.

end module arrays_store_useful