module theta_grid_params
  implicit none

  public :: init_theta_grid_params, init_trin_geo
  public :: wnml_theta_grid_params

  real, public :: rhoc, rmaj, r_geo, eps, epsl
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
         ntheta, nperiod, kp, asym, asympri, btor_slab

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
    in_file = input_unit_exist("theta_grid_parameters", exist)
    if (exist) read (unit=in_file, nml=theta_grid_parameters)

    if (kp > 0.) pk = 2.*kp

  end subroutine read_parameters

  subroutine wnml_theta_grid_params(unit)
   implicit none
   integer:: unit
       if (.not.exist) return
       write (unit, *)
       write (unit, fmt="(' &',a)") "theta_grid_parameters"
       write (unit, fmt="(' ntheta =  ',i4)") ntheta
       write (unit, fmt="(' nperiod = ',i4)") nperiod
       write (unit, fmt="(' rhoc =    ',f7.4)") rhoc
       write (unit, fmt="(' Rmaj =    ',f7.4)") rmaj
       write (unit, fmt="(' R_geo =   ',f7.4)") r_geo
       write (unit, fmt="(' eps =     ',f7.4)") eps
       write (unit, fmt="(' epsl =    ',f7.4)") epsl
       write (unit, fmt="(' qinp =    ',f7.4)") qinp
       write (unit, fmt="(' shat =    ',f7.4)") shat
       write (unit, fmt="(' alpmhd =  ',f7.4)") alpmhd
       write (unit, fmt="(' pk =      ',f7.4)") pk
       write (unit, fmt="(' kp =      ',f7.4)") kp
       write (unit, fmt="(' shift =   ',f7.4)") shift
       write (unit, fmt="(' akappa =  ',f7.4)") akappa
       write (unit, fmt="(' akappri = ',f7.4)") akappri
       write (unit, fmt="(' tri =     ',f7.4)") tri
       write (unit, fmt="(' tripri =  ',f7.4)") tripri
       write (unit, fmt="(' asym =     ',f7.4)") asym
       write (unit, fmt="(' asympri =  ',f7.4)") asympri
       write (unit, fmt="(' btor_slab =',f7.4)") btor_slab
       write (unit, fmt="(' /')")
  end subroutine wnml_theta_grid_params

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
