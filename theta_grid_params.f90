module theta_grid_params
  implicit none

  public :: init_theta_grid_params, init_trin_geo

  real, public :: rhoc, rmaj, r_geo, eps, epsl, drhotor2dr
  real, public :: qinp, shat, alpmhd, pk, shift, akappa, akappri, tri, tripri
  real, public :: asym, asympri, btor_slab
  real, public :: betaprim

  integer, public :: ntheta, nperiod

  private

  logical :: trin_flag = .false.
  real :: rhoc_trin, qval_trin, shat_trin, rgeo_trin, rmaj_trin, shift_trin
  real :: kappa_trin, kappri_trin, tri_trin, tripri_trin, betaprim_trin
  real :: kp = -1.
  logical :: exist

contains

  subroutine init_theta_grid_params
    implicit none
    logical, save :: initialized = .false.

    if (initialized) return
    initialized = .true.
    call read_parameters
    if (trin_flag) call reinit_theta_grid_params (rhoc_trin, qval_trin, &
         shat_trin, rgeo_trin, rmaj_trin, kappa_trin, kappri_trin, tri_trin, tripri_trin, &
         shift_trin, betaprim_trin)
  end subroutine init_theta_grid_params

  subroutine read_parameters
    use file_utils, only: input_unit, input_unit_exist
    implicit none
    integer :: in_file
!CMR,2/2/2011: add btor_slab
! btor_slab = btor/bpol defines direction of a flow relative to B in slab 
! geometry, where flow is by definition in the toroidal direction.
    namelist /theta_grid_parameters/ rhoc, rmaj, r_geo, eps, epsl, &
         qinp, shat, alpmhd, pk, shift, akappa, akappri, tri, tripri, &
         ntheta, nperiod, kp, asym, asympri, btor_slab, drhotor2dr

    rhoc = 0.5
    rmaj = 3.0
    r_geo = 3.0
    eps = 0.3
    epsl = 0.3
    qinp = 1.5
    shat = 0.75
    pk = 0.3
    shift = 0.0
    akappa = 1.0
    akappri = 0.0
    tri = 0.0
    tripri = 0.0
    asym = 0.0
    asympri = 0.0
    btor_slab = 0.0
    betaprim = 0.0
    ntheta = 24
    nperiod = 2
    drhotor2dr = 1.0

    in_file = input_unit_exist("theta_grid_parameters", exist)
    if (exist) read (unit=in_file, nml=theta_grid_parameters)

    if (kp > 0.) pk = 2.*kp

  end subroutine read_parameters

  subroutine reinit_theta_grid_params (rhoc_in, qval_in, shat_in, rgeo_in, rmaj_in, &
       kappa_in, kappri_in, tri_in, tripri_in, shift_in, betaprim_in)

    implicit none

    real, intent (in) :: rhoc_in, qval_in, shat_in, rgeo_in, rmaj_in, kappa_in, tri_in
    real, intent (in) :: kappri_in, tripri_in, shift_in, betaprim_in

    rhoc = rhoc_in
    qinp = qval_in
    shat = shat_in
    rmaj = rmaj_in
    r_geo = rgeo_in
    akappa = kappa_in
    akappri = kappri_in
    tri = tri_in
    tripri = tripri_in
    shift = shift_in
    betaprim = betaprim_in

  end subroutine reinit_theta_grid_params

  subroutine init_trin_geo (rhoc_in, qval_in, shat_in, rgeo_in, rmaj_in, &
       kappa_in, kappri_in, tri_in, tripri_in, shift_in, betaprim_in)

    implicit none

    real, intent (in) :: rhoc_in, qval_in, shat_in, rgeo_in, rmaj_in, kappa_in, tri_in
    real, intent (in) :: kappri_in, tripri_in, shift_in, betaprim_in

    trin_flag = .true.

    rhoc_trin = rhoc_in
    qval_trin = qval_in
    shat_trin = shat_in
    rgeo_trin = rgeo_in
    rmaj_trin = rmaj_in
    kappa_trin = kappa_in
    kappri_trin = kappri_in
    tri_trin = tri_in
    tripri_trin = tripri_in
    shift_trin = shift_in
    betaprim_trin = betaprim_in

  end subroutine init_trin_geo

end module theta_grid_params
