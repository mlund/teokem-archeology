! =============================================================================
! Common Blocks for Monte Carlo Simulation of Ionic Solutions
! All variables are explicitly typed (no reliance on implicit typing rules)
! =============================================================================

! =============================================================================
! Parameters
! =============================================================================
integer, parameter :: max_ions    = 600  ! Maximum number of ions in simulation
integer, parameter :: max_species = 10   ! Maximum number of ionic species


! =============================================================================
! Character Variables
! =============================================================================
! character(len=16) :: title  ! Simulation title (currently unused)


! =============================================================================
! Particle Position and Property Arrays
! =============================================================================
double precision :: particle_x(max_ions)         ! X coordinates of all particles
double precision :: particle_y(max_ions)         ! Y coordinates of all particles
double precision :: particle_z(max_ions)         ! Z coordinates of all particles
double precision :: particle_charge(max_ions)    ! Charge of each particle
double precision :: displacement_max(max_ions)   ! Maximum displacement for MC moves per particle


! =============================================================================
! Energy Arrays
! =============================================================================
double precision :: energy_matrix(max_ions, max_ions)  ! Pair interaction energy matrix
double precision :: trial_energy(max_ions)             ! Energy for trial particle configuration


! =============================================================================
! Species Properties Arrays
! =============================================================================
! Species properties matrix: (species, property_index)
! property_index: 1=species_id, 2=number_of_particles, 3=radius, 4=charge, 5=displacement
double precision :: species_properties(max_species, 5)

! Hard core interaction matrix (squared distances)
double precision :: hard_core_distance_sq(max_ions, max_species)

! Average concentration per species
double precision :: species_concentration(max_species)

! Species interaction matrix
double precision :: species_interaction_matrix(max_species, max_species)

! Contact correlation function g(r) at contact
double precision :: contact_correlation(25, max_species, max_species)

! Particle species assignment
integer :: particle_species(max_ions)  ! Species index for each particle


! =============================================================================
! Physical Constants (read-only parameters)
! =============================================================================
double precision, parameter :: pi = acos(-1.0d0)                 ! Mathematical constant π
double precision, parameter :: avogadro_number = 6.0223d23       ! Avogadro's number (N_A)
double precision, parameter :: gas_constant = 8.314d0            ! Gas constant R (kJ/mol/K)
double precision, parameter :: vacuum_permittivity = 8.85418782d-12  ! Vacuum permittivity (ε_0) in F/m
double precision, parameter :: elementary_charge = 1.60219d-19   ! Elementary charge (e) in C

! =============================================================================
! System Parameters (set from input or calculated)
! =============================================================================
double precision :: dielectric_constant     ! Relative dielectric constant (ε_r)
double precision :: temperature             ! Temperature (K)
double precision :: energy_conversion_factor ! Energy conversion factor for Coulomb interactions
double precision :: beta_inverse_temp       ! β = -1/(kT) for Metropolis criterion


! =============================================================================
! Trial Move Variables (for Monte Carlo moves)
! =============================================================================
double precision :: trial_x  ! Trial X position for particle move
double precision :: trial_y  ! Trial Y position for particle move
double precision :: trial_z  ! Trial Z position for particle move


! =============================================================================
! Simulation Box Parameters
! =============================================================================
double precision :: box_size          ! Simulation box size (cubic box)
double precision :: box_half          ! Half of box size (for periodic boundaries)
double precision :: box_half_inverse  ! Inverse of half box size (for efficiency)


! =============================================================================
! Energy Tracking Variables
! =============================================================================
double precision :: total_coulomb_energy  ! Total Coulomb energy of system
! double precision :: energy_check_parameter ! Parameter for energy consistency checks (currently unused)


! =============================================================================
! Pressure and Chemical Potential Variables
! =============================================================================
double precision :: pressure_widom_wall             ! Widom pressure at wall
double precision :: pressure_widom_midplane         ! Widom pressure at midplane
double precision :: pressure_widom_wall_variance    ! Variance of Widom wall pressure
double precision :: pressure_widom_midplane_variance ! Variance of Widom midplane pressure
double precision :: pressure_collision_average      ! Average collision pressure
double precision :: pressure_collision_variance     ! Variance of collision pressure


! =============================================================================
! Monte Carlo Loop Counters and Simulation Parameters
! =============================================================================
integer :: mc_steps_inner   ! Number of inner Monte Carlo steps
integer :: mc_steps_middle  ! Number of middle Monte Carlo steps
integer :: mc_steps_outer   ! Number of outer Monte Carlo steps (macro steps)
integer :: current_macro_step ! Current macro step counter

integer :: num_species      ! Total number of ionic species
integer :: current_species  ! Index of current species being processed
integer :: num_particles    ! Total number of particles in simulation

integer :: num_widom_insertions  ! Number of Widom test particle insertions
integer :: widom_interval        ! Interval between Widom calculations
integer :: measurement_location  ! Location for chemical potential measurement (0=bulk, 1=wall, 2=midplane)

integer :: current_particle  ! Index of current particle being moved
integer :: overlap_status    ! Overlap status flag (0=no overlap, non-zero=overlap)


! =============================================================================
! Configuration and I/O Parameters
! =============================================================================
integer :: initial_config_type  ! Initial configuration type (0=read from file, 1=random)
integer :: random_seed          ! Seed for random number generator


! =============================================================================
! File Unit Numbers
! =============================================================================
integer :: unit_conf    ! File unit for configuration file
integer :: unit_input   ! File unit for input file
integer :: unit_output  ! File unit for main output file
integer :: unit_misc    ! File unit for miscellaneous output
integer :: unit_macro   ! File unit for macro step output


! =============================================================================
! Common Block Declaration
! This preserves the memory layout for compatibility with legacy code
! NOTE: Physical constants (pi, avogadro_number, gas_constant, vacuum_permittivity,
!       elementary_charge) are now parameters and cannot be in common blocks
! =============================================================================
common /arr1/ energy_matrix, trial_energy, &
              species_concentration, species_interaction_matrix, &
              contact_correlation, &
              particle_x, particle_y, particle_z, particle_charge, displacement_max, &
              hard_core_distance_sq, species_properties, particle_species, &
              dielectric_constant, energy_conversion_factor, beta_inverse_temp, temperature, &
              trial_x, trial_y, trial_z, total_coulomb_energy, &
              pressure_widom_wall, pressure_widom_midplane, &
              pressure_widom_wall_variance, pressure_widom_midplane_variance, &
              pressure_collision_average, pressure_collision_variance, &
              box_size, box_half, box_half_inverse, &
              mc_steps_inner, mc_steps_middle, mc_steps_outer, current_macro_step, &
              num_species, current_species, num_particles, &
              num_widom_insertions, widom_interval, measurement_location, &
              current_particle, overlap_status, random_seed, initial_config_type, &
              unit_conf, unit_input, unit_output, unit_misc, unit_macro
