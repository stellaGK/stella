module parallel_dynamics

  use debug_flags, only: debug => parallel_dynamics_debug
  
  implicit none

  public :: init_parallel_dynamics
  public :: advance_parallel_dynamics
  public :: finish_parallel_dynamics

  private

  integer, dimension (:, :, :, :, :), allocatable :: departure_point_izext, departure_point_iv
  logical, dimension (:, :, :, :, :), allocatable :: departure_point_outside_grid
  
contains

  subroutine init_parallel_dynamics

    implicit none

    ! allocate arrays needed during time advance
    if (debug) write (*, *) 'parallel_dynamics::init_parallel_dynamics::allocate_arrays'
    call allocate_arrays

    ! for every grid point in the 2D (zext, vpa) domain, follow the corresponding
    ! lowest-order characteristic -- which has constant particle kinetic energy --
    ! back in time a distance dt to find the 'departure point' in the (zext, vpa) space
    if (debug) write (*, *) 'parallel_dynamics::init_parallel_dynamics::find_departure_points'
    call find_departure_points
    
  end subroutine init_parallel_dynamics

  subroutine allocate_arrays

    use parameters_kxky_grids, only: nakx
    use zgrid, only: nzgrid, ntubes
    use vpamu_grids, only: nvpa
    use stella_layouts, only: kymus_lo
    
    implicit none

    if (.not. allocated(departure_point_izext)) allocate(departure_point_izext(nakx, -nzgrid:nzgrid, ntubes, nvpa, kymus_lo%llim_proc:kymus_lo%ulim_alloc))
    if (.not. allocated(departure_point_iv)) allocate(departure_point_iv(nakx, -nzgrid:nzgrid, ntubes, nvpa, kymus_lo%llim_proc:kymus_lo%ulim_alloc))
    if (.not. allocated(departure_point_outside_grid)) allocate(departure_point_outside_grid(nakx, -nzgrid:nzgrid, ntubes, nvpa, kymus_lo%llim_proc:kymus_lo%ulim_alloc))
    
  end subroutine allocate_arrays

  ! find_departure_points uses the fact that the particle kinetic energy is an
  ! approximate constant of the motion to follow a particle's characteristic
  ! backward in time an amount delt to find the particle's 'departure point'
  ! in the (zext, vpa) plane.
  ! In particular, dz/dt = vpa * (bhat . grad z) = +/- (bhat . grad z) * sqrt(2/m) * sqrt(E - mu * B(z)),
  ! with E and mu constants of the motion.
  ! Thus, the time to cross one grid cell of width dz centred at z_i is dt_i = dz / (bhat . grad z) / sqrt(2/m) / sqrt(E - mu * B(z_i)).
  ! This subroutine repeats this process until delt = sum_i dt_i, giving the departure z location.
  ! Once the departure z is known, the departure vpa is sqrt(2/m) * sqrt(E - mu * B(z_departure))
  subroutine find_departure_points

    use zgrid, only: nzgrid, ntubes
    use vpamu_grids, only: nvpa, vpa, mu
    use stella_layouts, only: kymus_lo
    use stella_layouts, only: imu_idx, iky_idx
    use geometry, only: bmag
    use extended_zgrid, only: nsegments, nzed_segment, neigen
    use extended_zgrid, only: map_to_iz_ikx_from_izext
    ! TMP FOR TESTING -- MAB
    use stella_time, only: code_dt
    use zgrid, only: delzed
    use geometry, only: b_dot_grad_z
    use vpamu_grids, only: nmu
    
    implicit none

    integer :: ikymus, iky, imu
    integer :: ikx, iz, it, iv
    integer :: izext, nz_ext, ie
    integer :: ia = 1  ! flux annulus is not supported
    integer, dimension (:), allocatable :: iz_from_izext, ikx_from_izext
    real, dimension (:, :), allocatable :: energy
    real, dimension (:), allocatable :: bmag_cell_centres
    
    allocate (energy(-nzgrid:nzgrid, nvpa))
    allocate (bmag_cell_centres(-nzgrid:nzgrid-1))
    
    ! obtain the normalised B-field at cell-centres, needed for the characteristic tracking algorithm
    do iz = -nzgrid, nzgrid-1
       bmag_cell_centres(iz) = 0.5 * (bmag(ia, iz) + bmag(ia, iz+1))
    end do
    
    ! the characteristic at a given (zext, vpa) value is determined by the particle's energy and mu
    ! energy = v^2 / vths^2 = energy(iv, imu, iz)
    ! the departure point will thus be a function of iv, imu, and izext (i.e., iz, ikx and it)
    do ikymus = kymus_lo%llim_proc, kymus_lo%ulim_proc
       imu = imu_idx(kymus_lo, ikymus)
       iky = iky_idx(kymus_lo, ikymus)
       ! calculate the particle kinetic energy as a function of zed and vpa at the given mu
       do iv = 1, nvpa
          do iz = -nzgrid, nzgrid
             energy(iz, iv) = vpa(iv)**2 + 2.0 * mu(imu) * bmag(ia, iz)
          end do
       end do
       do it = 1, ntubes
          do ie = 1, neigen(iky)
             ! nz_ext is the number of grid points in the extended zed domain
             nz_ext = nsegments(ie, iky) * nzed_segment + 1
             ! obtain the mapping to iz and ikx from the extended zed domain
             allocate (iz_from_izext(nz_ext)) ; allocate(ikx_from_izext(nz_ext))
             call map_to_iz_ikx_from_izext(iky, ie, iz_from_izext, ikx_from_izext)
             do iv = 1, nvpa
                do izext = 1, nz_ext
                   ! get the ikx and iz corresponding to this izext
                   ikx = ikx_from_izext(izext)
                   iz = iz_from_izext(izext)
                   call calculate_zed_vpa_departure_idx(sign(1.0, vpa(iv)), energy(iz_from_izext(izext), iv), mu(imu), &
                        izext, nz_ext, iz_from_izext, ikx_from_izext, departure_point_izext(ikx, iz, it, iv, ikymus), &
                        departure_point_iv(ikx, iz, it, iv, ikymus), departure_point_outside_grid(ikx, iz, it, iv, ikymus), iv)
                   if (imu == nmu/2) write (101,*) izext, iv, departure_point_izext(ikx, iz, it, iv, ikymus) - izext, departure_point_iv(ikx, iz, it, iv, ikymus) - iv, departure_point_outside_grid(ikx, iz, it, iv, ikymus)
!                   if (iv == 4 .and. imu == nmu .and. izext == nz_ext/2) write (*,*) 'departure_point_izext: ', izext, departure_point_izext(ikx, iz, it, iv, ikymus), departure_point_iv(ikx, iz, it, iv, ikymus), vpa(iv) * b_dot_grad_z(ia, iz) * code_dt / delzed(0)
                end do
                if (imu == nmu/2) write (101,*)
             end do
             deallocate(iz_from_izext, ikx_from_izext)
          end do
       end do
             
    end do

    deallocate (energy)
    deallocate (bmag_cell_centres)

  end subroutine find_departure_points

  subroutine calculate_zed_vpa_departure_idx(sign_vpa_initial, energy, mu, initial_izext, nz_ext, &
       iz_from_izext, ikx_from_izext, departure_izext, departure_iv, departure_outside_grid, iv)

    use stella_time, only: code_dt
    use zgrid, only: delzed
    use vpamu_grids, only: vpa, nvpa
    use geometry, only: bmag
    
    implicit none

    real, intent (in) :: sign_vpa_initial, energy, mu
    integer, intent (in) :: initial_izext, nz_ext
    integer, dimension (:), intent (in) :: iz_from_izext, ikx_from_izext
    integer, intent (out) :: departure_izext, departure_iv
    logical, intent (out) :: departure_outside_grid
    ! TMP FOR TESTING -- MAB
    integer, intent (in) :: iv
    
    integer :: iz, izext, izext_next
    integer :: ia = 1  ! flux annulus not supported
    integer :: sign_vpa
    real :: dt_cumulative, dt
    real :: dz_bounce, dt_bounce
    real :: dz_midpoint, dt_midpoint
    real :: vpa_departure
    logical :: bounce_point_in_cell
    
    ! dt_cumulative is the total time elapsed while following the characteristic
    ! across z-grid cells; initialise to zero.
    dt_cumulative = 0.0
    ! departure_outside_grid will be set to .true. if the characteristic leaves the zed domain
    departure_outside_grid = .false.
    
    ! the initial sign of the parallel velocity has been passed in as an argument;
    ! this will tell us in which direction the characteristic starts moving
    sign_vpa = int(sign_vpa_initial)
    
    ! set the initial izext value to the input value
    izext = initial_izext
    ! izext_next is the izext index towards which the characteristic is moving as we
    ! track it backwards in time
    izext_next = izext - sign_vpa
    
    ! track the characteristic across z-grid cells until the accumulated time taken
    ! for the particle trajectory exceeds the time step size or leaves the zed domain
    do while (dt_cumulative < code_dt .and. .not. departure_outside_grid)
       ! determine if the characteristic is due to leave the zed domain in this iteration
       if (izext_next < 1 .or. izext_next > nz_ext) then
          ! if the characteristic is leaving the zed domain, no need to keep tracking it so exit the do loop
          departure_outside_grid = .true. ; exit
       end if
       ! determine if the characteristic will change direction before crossing the next grid cell;
       ! i.e., is there a bounce point in this grid cell
       call check_for_bounce_point(izext_next, iz_from_izext, energy, mu, bounce_point_in_cell)
       ! if there is a bounce point in the cell, calculate the time needed to bounce and return
       ! to the current izext location
       if (bounce_point_in_cell) then
          ! find the location in zed of the bounce point
          call calculate_distance_to_bounce_point(izext, izext_next, iz_from_izext, sign_vpa, energy, mu, dz_bounce)
          ! find the time it takes to reach the bounce point
          call calculate_time_to_move_dz(izext, izext_next, iz_from_izext, sign_vpa, energy, mu, dz_bounce, dt_bounce)
          ! the time taken to return to the starting izext is double the time taken to reach the bounce point (symmetry)
          dt = 2 * dt_bounce
          ! check to see if adding dt to the cumulative time taken along the characteristic
          ! makes the cumulative time larger than the time step size
          if (dt + dt_cumulative > code_dt) then
             ! dz_midpoint is the distance to the centre of the cell in zed
             dz_midpoint = 0.5 * delzed(iz_from_izext(izext))
             ! if the bounce point is closer to z(izext) than z(izext_next), then the
             ! departure point must also be closer to z(izext) than z(izext_next)
             if (dz_bounce < dz_midpoint) then
                departure_izext = izext
             else if (dt_bounce + dt_cumulative > code_dt) then
                ! if code_dt is passed before reaching the bounce point, particle
                ! does not bounce; check to see if it reaches the cell centre
                ! before code_dt is reached

                ! use nearest-neighbour interpolation to place departure point on a zed grid point
                call find_nearest_zgrid_index(izext, izext_next, iz_from_izext, sign_vpa, energy, mu, dt_cumulative, departure_izext)
             else
                ! if the bounce point is beyond the cell centre,
                ! and the particle has time to bounce before code_dt is reached,
                ! then the time needed to bounce and make it back to the cell centre is
                ! dt - dt_midpoint; if dt_cumulative + dt - dt_midpoint < code_dt, then nearest grid index is izext;
                ! else, the nearest index is izext_next

                ! find dt_midpoint = the time taken to reach the mid-point of the cell
                call calculate_time_to_move_dz(izext, izext_next, iz_from_izext, sign_vpa, energy, mu, dz_midpoint, dt_midpoint)
                ! use nearest-neighbour interpolation to place departure point on a zed grid point
                if (dt_cumulative + dt - dt_midpoint < code_dt) then
                   departure_izext = izext
                else
                   departure_izext = izext_next
                end if
             end if
             ! TMP FOR TESTING -- MAB
             ! this should possibly be made a bit better -- ensures that no info beyond bounce point is used
             departure_izext = izext
             ! check to see if the particle has bounced before code_dt is reached;
             ! if so, the sign of the particle parallel velocity changes
             if (dt_bounce + dt_cumulative < code_dt) sign_vpa = -sign_vpa
          else
             ! after bouncing, the sign of the particle parallel velocity changes
             sign_vpa = -sign_vpa
          end if
          ! note that the updated izext is unchanged (only the sign of the velocity has changed)
       else
          ! if there is no bounce point, calculate the time it takes to cross the cell
          call calculate_time_to_move_dz(izext, izext_next, iz_from_izext, sign_vpa, energy, mu, delzed(iz_from_izext(izext)), dt)
          ! check to see if adding dt to the cumulative time taken along the characteristic
          ! makes the cumulative time larger than the time step size
          if (dt + dt_cumulative > code_dt) then
             ! use nearest-neighbour interpolation to place departure point on a zed grid point
             call find_nearest_zgrid_index(izext, izext_next, iz_from_izext, sign_vpa, energy, mu, dt_cumulative, departure_izext)
          else
             ! update the izext of the particle as it moves (backward) along the characteristic
             izext = izext_next
          end if
       end if
          
       ! update the cumulative time taken to track the characteristic
       dt_cumulative = dt_cumulative + dt

       if (dt_cumulative > code_dt) then
          ! obtain the corresponding vpa value for the departure point using E, mu constant
          vpa_departure = sign_vpa * sqrt(energy - 2.0 * mu * bmag(ia, iz_from_izext(departure_izext)))
          if (abs(vpa_departure) > vpa(nvpa)) departure_outside_grid = .true.
          ! obtain the associated iv index
          departure_iv = minloc(abs(vpa - vpa_departure), dim=1)
!          if (iv == 49) then
!             write (*,*) 'vpa_departure', iv, vpa_departure, sign_vpa * sqrt(energy - 2.0 * mu * bmag(ia, iz_from_izext(izext))), departure_izext - izext, izext, departure_outside_grid
!          end if
       end if
          
       ! izext_next is the izext index towards which the characteristic is moving as we track it backwards in time
       izext_next = izext - sign_vpa
    end do
      
  end subroutine calculate_zed_vpa_departure_idx

  subroutine check_for_bounce_point(izext_next, iz_from_izext, energy, mu, bounce_point_in_cell)

    use geometry, only: bmag
    
    implicit none

    integer, intent (in) :: izext_next
    integer, dimension (:), intent (in) :: iz_from_izext
    real, intent (in) :: energy, mu
    logical, intent (out) :: bounce_point_in_cell
    
    integer :: iz_next
    integer :: ia = 1  ! flux annulus is not supported
    real :: parallel_kinetic_energy
    
    ! initialise bounce_point_in_cell to .false. and change to .true.
    ! if a bounce point is found in the cell
    bounce_point_in_cell = .false.

    ! iz_next is the iz index towards which the characteristic is moving as we track it backwards
    ! in time
    iz_next = iz_from_izext(izext_next)
    ! parallel_kinetic_energy is the contribution to the particle kinetic energy due to
    ! parallel motion; it is evaluated at iz_next
    parallel_kinetic_energy = energy - 2.0 * mu * bmag(ia, iz_next)
    ! as parallel_kinetic_energy is positive at iz, a change in sign indicates the
    ! presence of a bounce point in the cell
    if (parallel_kinetic_energy < 0.0) bounce_point_in_cell = .true.
    
  end subroutine check_for_bounce_point

  subroutine calculate_distance_to_bounce_point(izext, izext_next, iz_from_izext, sign_vpa, energy, mu, dz)

    use geometry, only: bmag
    use zgrid, only: delzed
    
    implicit none

    integer, intent (in) :: izext, izext_next, sign_vpa
    integer, dimension (:), intent (in) :: iz_from_izext
    real, intent (in) :: energy, mu
    real, intent (out) :: dz

    integer :: iz, iz_next
    integer :: ia = 1  ! flux annulus is not supported
    real :: dB_dz, parallel_kinetic_energy

    ! this is the starting iz index for this iteration of the particle tracking
    iz = iz_from_izext(izext)
    ! this is the iz index towards which the characteristic is moving as we track it backward in time
    iz_next = iz_from_izext(izext_next)
    ! this is the contribution to the kinetic energy at the starting iz index coming from parallel motion
    parallel_kinetic_energy = energy - 2.0 * mu * bmag(ia, iz)

    ! bounce occurs when E = 2 * mu * B(z)
    ! using B(z) = B(z(izext)) + (z - z(izext)) * (dB_dz)_izext,
    ! with (dB_dz)_izext = (B(z(izext)) - B(z(izext_next))) / (z(izext) - z(izext_next)),
    ! gives E - 2 * mu * B(z(izext)) = (z - z(izext)) * (dB_dz)_izext.
    ! dz = | z - z(izext) |
    ! only need absolute value of dB/dz to get dz (which is defined to be positive)
    dB_dz = abs(bmag(ia, iz) - bmag(ia, iz_next)) / delzed(iz)
    dz = parallel_kinetic_energy / dB_dz
    
  end subroutine calculate_distance_to_bounce_point
  
  ! calculate the time taken for a particle to move the distance dz
  ! along the characteristic defined by energy = particle kinetic energy
  ! and mu = particle magnetic moment
  subroutine calculate_time_to_move_dz(izext, izext_next, iz_from_izext, sign_vpa, energy, mu, dz, dt)

    use geometry, only: bmag, b_dot_grad_z
    use zgrid, only: delzed
    
    implicit none

    integer, intent (in) :: izext, izext_next, sign_vpa
    integer, dimension (:), intent (in) :: iz_from_izext
    real, intent (in) :: energy, mu, dz
    real, intent (out) :: dt

    integer :: iz, iz_next
    integer :: ia = 1  ! flux annulus is not supported
    real :: parallel_kinetic_energy, dB_dz

    ! this is the starting iz index for this iteration of the particle tracking
    iz = iz_from_izext(izext)
    ! this is the iz index towards which the characteristic is moving as we track it backward in time
    iz_next = iz_from_izext(izext_next)
    ! this is the contribution to the kinetic energy at the starting iz index coming from parallel motion
    parallel_kinetic_energy = energy - 2.0 * mu * bmag(ia, iz)
    ! we use linear interpolation to fill in B(z) between grid points:
    ! B(z) = B(z(izext)) + (z - z(izext)) * (dB/dz)_izext,
    ! with (dB/dz)_izext = ( B(z(izext)) - B(z(izext_next)) ) / ( z(izext) - z(izext_next) ).
    ! sign_vpa necessary because delzed = | z(izext) - z(izext_next) |
    dB_dz = sign_vpa * (bmag(ia, iz) - bmag(ia, iz_next)) / delzed(iz)
    ! integrate dz/dt = sign_vpa * (b . grad z) * sqrt(E - 2 * mu * B(z))
    ! to get the time dt taken to move a distance dz
    dt = 2.0 * abs(( sqrt(parallel_kinetic_energy + mu * sign_vpa * dz * dB_dz) - sqrt(parallel_kinetic_energy) )) &
         / abs(mu * b_dot_grad_z(ia, iz) * dB_dz)

  end subroutine calculate_time_to_move_dz

  subroutine find_nearest_zgrid_index(izext, izext_next, iz_from_izext, sign_vpa, energy, mu, dt_cumulative, departure_izext)

    use zgrid, only: delzed
    use stella_time, only: code_dt
    
    implicit none

    integer, intent (in) :: izext, izext_next, sign_vpa
    integer, dimension (:), intent (in) :: iz_from_izext
    real, intent (in) :: energy, mu, dt_cumulative
    integer, intent (out) :: departure_izext
    
    real :: dz_midpoint, dt_midpoint
    
    ! dz_midpoint is the distance to the centre of the cell in zed
    dz_midpoint = 0.5 * delzed(iz_from_izext(izext))
    
    ! find dt_midpoint = the time taken to reach the mid-point of the cell
    call calculate_time_to_move_dz(izext, izext_next, iz_from_izext, sign_vpa, energy, mu, dz_midpoint, dt_midpoint)
    
    ! use nearest-neighbour interpolation to place departure point on a zed grid point
    if (dt_midpoint + dt_cumulative > code_dt) then
       departure_izext = izext
    else
       departure_izext = izext_next
    end if
    
  end subroutine find_nearest_zgrid_index
  
  subroutine advance_parallel_dynamics (pdf, phi)

    use finite_differences, only: second_order_centered
    use zgrid, only: nzgrid, ntubes
    use zgrid, only: delzed
    use vpamu_grids, only: nvpa, vpa
    use vpamu_grids, only: maxwell_vpa, maxwell_mu
    use species, only: spec
    use stella_time, only: code_dt
    use stella_layouts, only: vmu_lo, kymus_lo
    use stella_layouts, only: iky_idx
    use redistribute, only: scatter, gather
    use dist_redistribute, only: kymus2vmus
    use arrays_dist_fn, only: g_kymus
    use extended_zgrid, only: nsegments, nzed_segment, neigen
    use extended_zgrid, only: map_to_extended_zgrid, map_from_extended_zgrid
    use extended_zgrid, only: map_to_iz_ikx_from_izext
    
    implicit none

    complex, dimension (:, :, -nzgrid:, :, vmu_lo%llim_proc:), intent (in out) :: pdf
    complex, dimension (:, :, -nzgrid:, :), intent (in out) :: phi

    integer :: ikymus, iky, imu, is
    integer :: ikx, iz, it, iv, ie
    integer :: izext, nz_ext
    integer :: ulim   ! dummy variable
    integer :: ia = 1 ! does not support flux annulus
    integer, dimension (:), allocatable :: iz_from_izext, ikx_from_izext
    complex, dimension (:, :), allocatable :: g_ext_1, g_ext_2
    complex, dimension (:), allocatable :: phi_ext, dphi_dz
    
    ! the input pdf is in the vmu_lo; re-map to work in kymus_lo
    call scatter(kymus2vmus, pdf, g_kymus)

    ! update g = <delta f> to account for motion along the lowest-order characteristic,
    ! for which the particle kinetic energy is constant
    do ikymus = kymus_lo%llim_proc, kymus_lo%ulim_proc
       ! map to the extended zed domain to ease calculations
       iky = iky_idx(kymus_lo, ikymus)
       do it = 1, ntubes
          do ie = 1, neigen(iky)
             ! nz_ext is the number of grid points in the extended zed domain
             nz_ext = nsegments(ie, iky) * nzed_segment + 1
             allocate (iz_from_izext(nz_ext))
             allocate (ikx_from_izext(nz_ext))
             ! g_ext_1 and g_ext_2 will contain slices of the pdf on the extended zed domain
             allocate (g_ext_1(nz_ext, nvpa)) ; allocate (g_ext_2(nz_ext, nvpa))
             do iv = 1, nvpa
                ! map from (kx, z, tube) to the extended zed domain; NB: ulim is a dummy argument
                call map_to_extended_zgrid(it, ie, iky, g_kymus(:, :, :, iv, ikymus), g_ext_1(:, iv), ulim)
             end do
             call map_to_iz_ikx_from_izext(iky, ie, iz_from_izext, ikx_from_izext)
             ! set the pdf(t+dt) at each (zext, vpa) grid point to be equal to
             ! the value of the pdf(t) at the (zext, vpa) that connects to it along the
             ! characteristic with constant particle kinetic energy
             ! NB: currently using a crude nearest-neighbour interpolation for departure point
             do iv = 1, nvpa
                do izext = 1, nz_ext
                   iz = iz_from_izext(izext)
                   ikx = ikx_from_izext(izext)
                   if (departure_point_outside_grid(ikx, iz, it, iv, ikymus)) then
                      g_ext_2(izext, iv) = 0.
                   else
                      g_ext_2(izext, iv) = g_ext_1(departure_point_izext(ikx, iz, it, iv, ikymus), departure_point_iv(ikx, iz, it, iv, ikymus))
                   end if
                end do
             end do
             ! phi_ext will contain a slice of the electrostatic potential on the extended zed domain
             allocate (phi_ext(nz_ext))
             ! dphi_dz will contain a slice of dphi/dz on the extended zed domain
             allocate (dphi_dz(nz_ext))
             ! map from (kx, z, tube) to the extended zed domain; NB: ulim is a dummy argument
             call map_to_extended_zgrid(it, ie, iky, phi(iky, :, :, :), phi_ext, ulim)
             ! compute dphi/dz using centered differences, with zero BCs
             call second_order_centered(1, phi_ext, delzed(0), dphi_dz)
             ! calculate the change in g over time dt due to the parallel acceleration of
             ! particles in the background Maxwellian population (by the parallel electric field);
             ! for now, use the parallel electric field at the departure point to obtain
             ! the acceleration
             do iv = 1, nvpa
                do izext = 1, nz_ext
                   iz = iz_from_izext(izext)
                   ikx = ikx_from_izext(izext)
                   g_ext_2(izext, iv) = g_ext_2(izext, iv) + code_dt * spec(is)%zstm &
                        * maxwell_vpa(iv, is) * maxwell_mu(ia, iz, imu, is) &
                        * vpa(departure_point_iv(ikx, iz, it, iv, ikymus)) * dphi_dz(departure_point_izext(ikx, iz, it, iv, ikymus))
                end do
             end do
             do iv = 1, nvpa
                ! map from the extended zed domain to (kx, z, tube); NB: ulim is a dummy argument
                call map_from_extended_zgrid(it, ie, iky, g_ext_2(:, iv), g_kymus(:, :, :, iv, ikymus))
             end do
             ! deallocate g_ext so that it can be re-allocated with different size
             deallocate (g_ext_1, g_ext_2, phi_ext, dphi_dz)
             deallocate (iz_from_izext, ikx_from_izext)
          end do
       end do
    end do
    
    ! the output pdf should be in the vmu_lo; re-map from kymus_lo to vmu_lo
    call gather(kymus2vmus, g_kymus, pdf)
    
  end subroutine advance_parallel_dynamics

  subroutine finish_parallel_dynamics

    implicit none

    if (allocated(departure_point_izext)) deallocate(departure_point_izext)
    if (allocated(departure_point_iv)) deallocate(departure_point_iv)
    if (allocated(departure_point_outside_grid)) deallocate(departure_point_outside_grid)
    
  end subroutine finish_parallel_dynamics
  
end module parallel_dynamics
