! =============================================================================
! =============================================================================
!
! MONTE CARLO SIMULATION OF MULTICOMPONENT ISOTROPIC IONIC SOLUTION
!
! =============================================================================
! =============================================================================

program bulk

   implicit none
   double precision :: random_value
   include 'bulk_f90.inc'

   ! Explicit interfaces for subroutines
   interface
      subroutine calculate_statistics(data_array, standard_deviation, mean_value, array_size)
         integer :: array_size
         double precision :: data_array(25), standard_deviation, mean_value
      end subroutine calculate_statistics

     subroutine calculate_statistics_per_species(data_matrix, standard_deviations, mean_values, num_samples, num_species_to_process)
         integer, parameter :: max_species = 10
         integer :: num_samples, num_species_to_process
         double precision :: data_matrix(25, max_species), mean_values(max_species), standard_deviations(max_species)
      end subroutine calculate_statistics_per_species

      subroutine initialize_random_configuration()
      end subroutine initialize_random_configuration

      subroutine evaluate_trial_move()
      end subroutine evaluate_trial_move

      subroutine recalculate_total_energy()
      end subroutine recalculate_total_energy

      subroutine calculate_collision_pressure()
      end subroutine calculate_collision_pressure

      subroutine calculate_collision_pressure1()
      end subroutine calculate_collision_pressure1

      subroutine calculate_collision_pressure2()
      end subroutine calculate_collision_pressure2

      subroutine calculate_collision_pressure3()
      end subroutine calculate_collision_pressure3

      subroutine calculate_widom_insertion()
      end subroutine calculate_widom_insertion

      subroutine calculate_widom_insertion1()
      end subroutine calculate_widom_insertion1

      subroutine calculate_widom_insertion2()
      end subroutine calculate_widom_insertion2

      subroutine calculate_widom_insertion3()
      end subroutine calculate_widom_insertion3

      subroutine initialize_rdf()
      end subroutine initialize_rdf

      subroutine accumulate_rdf()
      end subroutine accumulate_rdf

      subroutine finalize_and_write_rdf()
      end subroutine finalize_and_write_rdf
   end interface

   ! Note: Variable names have been updated to be more descriptive
   ! See bulk_f90.inc for the complete variable declarations

   ! Local variables for Monte Carlo simulation
   integer :: step_inner, step_middle                        ! Loop counters for inner and middle MC loops
   integer :: total_hardcore_rejections                      ! Total hard core overlap rejections
   integer :: total_energy_rejections                        ! Total energy-based rejections
   integer :: total_accepted_moves                           ! Total accepted Monte Carlo moves
   integer :: total_mc_steps                                 ! Total Monte Carlo steps attempted
   integer :: species_index, particle_index                  ! Loop indices for species and particles
   integer :: loop_i, loop_j, loop_k, loop_l                 ! Generic loop counters
   integer :: particle_count_temp                            ! Temporary particle count
   integer :: total_configurations                           ! Total number of configurations
   integer :: total_widom_tests                              ! Total Widom test particle insertions
   integer :: moves_per_species_total(max_species)           ! Total moves attempted per species
   integer :: moves_per_species_accepted(max_species)        ! Accepted moves per species
   integer :: moves_per_species_energy_rejected(max_species) ! Energy rejected moves per species
   integer :: moves_per_species_hardcore_rejected(max_species) ! Hard core rejected moves per species

   double precision :: energy_accumulator_current(2)         ! Energy accumulator for current macro step
   double precision :: energy_accumulator_total(2)           ! Total energy accumulator across all macro steps
   double precision :: energy_per_macrostep(25)              ! Energy stored for each macro step
   double precision :: time_array(10)                        ! Time measurements array (currently unused)
   double precision :: energy_difference                     ! Energy difference for MC acceptance
   double precision :: boltzmann_factor                      ! Boltzmann factor for MC acceptance
   double precision :: energy_variance                       ! Variance of energy measurements
   double precision :: temp_kelvin                           ! Temperature in Kelvin
   double precision :: pressure_excess_energy                ! Excess pressure from energy
   double precision :: pressure_excess_variance              ! Variance of excess pressure
   double precision :: pressure_ideal                        ! Ideal gas pressure contribution
   double precision :: pressure_total                        ! Total system pressure
   double precision :: pressure_total_variance               ! Variance of total pressure
   double precision :: energy_consistency_check              ! Relative energy consistency check
   double precision :: macro_step_inverse                    ! Inverse of number of macro steps (1/N)
   double precision :: initial_energy_scaled                 ! Initial configuration energy (scaled to kJ/mol)

   data energy_accumulator_total/2*0.0/
   data energy_per_macrostep/25*0.0/
   data moves_per_species_total, moves_per_species_accepted, &
      moves_per_species_energy_rejected, moves_per_species_hardcore_rejected/40*0/
   data total_hardcore_rejections, total_energy_rejections, &
      total_accepted_moves, total_mc_steps/4*0/
   data time_array/10*0.0/

   unit_macro = 12
   unit_conf = 13
   unit_input = 15
   unit_output = 6

   write (unit_output, '(/, "************************************************************", &
   &"************************", /)')
   write (unit_output, '("MONTE CARLO SIMULATION OF MULTICOMPONENT ISOTROPIC    ", &
   &"IONIC SYSTEM  ", /, /, "Division of Theoretical Chemistry, Lund University", /)')
   write (unit_output, '(/, "************************************************************", &
   &"************************", /)')

   ! ==========================================================================
   ! Physical constants are now defined as parameters in bulk_f90.inc
   ! ==========================================================================

   ! ==========================================================================
   ! Sample units are angstrom and elementary charge
   ! Kilojoules are used for a mole
   !
   ! MEANING OF THE INPUT VARIABLES
   !
   ! initial_config_type        = 0   from file
   ! initial_config_type        = 1   from slump
   ! num_species      = number of ionic components
   ! hion        = array contains the input info on the ions
   !               (1=nr, 2=number, 3=radius, 4=charge, 5=displacement)
   ! ny(1,2,3)   = number of MC loops
   ! temperature      = temperature
   ! dielectric_constant        = relative dielectric constant
   ! box_size        = box_sizesize
   ! widom_interval      = a widom/collision calc is performed every widom_intervalcycle
   ! num_widom_insertions      = number of widom particles inserted each calc.
   ! measurement_location       = where is the chem. pot. measured
   !               0 = averaged over the cell
   !
   ! ==========================================================================
   ! ==========================================================================

   open (unit_conf, file='bulk.conf', form='formatted')
   open (unit_input, file='bulk.inp', form='formatted')
   open (unit_macro, file='bulk.macro', form='formatted')

   read (unit_input, *)
   read (unit_input, *) num_species
   read (unit_input, *)
   read (unit_input, *)
   read(unit_input, *) (species_properties(loop_l, 1), species_properties(loop_l, 2), species_properties(loop_l, 3), species_properties(loop_l, 4), species_properties(loop_l, 5), loop_l = 1, num_species)
   read (unit_input, *)
   read (unit_input, *)
   read (unit_input, *) mc_steps_inner, mc_steps_middle, mc_steps_outer
   read (unit_input, *)
   read (unit_input, *)
   read (unit_input, *) initial_config_type, random_seed
   read (unit_input, *)
   read (unit_input, *)
   read (unit_input, *) box_size
   read (unit_input, *)
   read (unit_input, *)
   read (unit_input, *) temperature, dielectric_constant
   read (unit_input, *)
   read (unit_input, *)
   read (unit_input, *) num_widom_insertions, widom_interval, measurement_location

   if (measurement_location > 2) measurement_location = 2
   if (random_seed == 0) random_seed = 7
   if (widom_interval == 0) widom_interval = 2*mc_steps_inner

   ! Initialize built-in random number generator with seed
   ! Note: The common block contains a variable named random_seed which conflicts
   ! with the intrinsic random_seed subroutine, so we use a helper subroutine
   call initialize_rng(random_seed)

   ! ==========================================================================
   ! Boxsize variables
   ! ==========================================================================

   box_half = box_size/2.0
   box_half_inverse = 1.0/box_half

   ! ==========================================================================
   ! Set up hard core and charge vectors
   !
   ! hc2v contains in the first column the particle type given
   ! by species_properties(i,1), whereas in the following columns the squared
   ! hard core distances between the particles are given.
   ! chv gives for every particle its charge
   ! num_particlesis the total number of particles
   ! ==========================================================================

   num_particles = 0

   do species_index = 1, num_species
      particle_count_temp = int(species_properties(species_index, 2))

      do particle_index = 1, particle_count_temp
         num_particles = num_particles + 1
         particle_species(num_particles) = species_index

         do loop_k = 1, num_species
            hard_core_distance_sq(num_particles, loop_k) = (species_properties(species_index, 3) + species_properties(loop_k, 3))**2
         end do

         particle_charge(num_particles) = species_properties(species_index, 4)
         displacement_max(num_particles) = species_properties(species_index, 5)
      end do
   end do

   do species_index = 1, num_species
      species_concentration(species_index) = (species_properties(species_index, 2))/(box_size**3)
   end do

   do particle_index = 1, num_particles
      trial_energy(particle_index) = 0.0
   end do

   total_configurations = mc_steps_inner*mc_steps_middle*mc_steps_outer

   if (mc_steps_outer > 25) mc_steps_outer = 25

   temp_kelvin = temperature
   beta_inverse_temp = -1.0/(gas_constant*temperature)

   ! ==========================================================================
   ! Read configuration from file
   ! ==========================================================================

   if (initial_config_type == 0) then
      read (unit_conf, *) (particle_x(loop_l), particle_y(loop_l), particle_z(loop_l), loop_l=1, num_particles)
      write (unit_output, *) 'Configuration read from file'
   end if

   if (initial_config_type == 1) call initialize_random_configuration
   call initialize_rdf
   call calculate_collision_pressure
   call calculate_widom_insertion

   ! ==========================================================================
   ! Write input data
   ! ==========================================================================

   total_widom_tests = num_widom_insertions*(int(mc_steps_inner/widom_interval)*mc_steps_middle*mc_steps_outer)

   write (unit_output, '(/, "VARIABLES FROM INPUT FILE", /, &
   &/, " number of particles        =", i10, /, &
   &/, " species   number   radius    charge   dp   average conc", /)') num_particles
   write (unit_output, '(i6, i10, f10.1, f10.1, f10.2, f10.6)') &
      (loop_l, int(species_properties(loop_l, 2)), species_properties(loop_l, 3), species_properties(loop_l, 4), &
       species_properties(loop_l, 5), species_concentration(loop_l)*1.0d27/avogadro_number, loop_l=1, num_species)
   write (unit_output, '(/, " dielectric constant         =", f10.1, /, &
   &" temperature                 =", f10.1, /, &
   &" box_sizesize (AA)               =", f10.1, /, &
   &" number of configurations    =", i6, " *", i6, " *", i6, " =", i8, /, &
   &" number of conf. per part    =", f10.1, /, &
   &" widom test per species      =", i10, /)') &
      dielectric_constant, temperature, box_size, mc_steps_inner, mc_steps_middle, mc_steps_outer, total_configurations, &
      dble(total_configurations/num_particles), total_widom_tests

   energy_conversion_factor = elementary_charge**2*1.0d10*avogadro_number/(dielectric_constant*vacuum_permittivity*4.0*pi)

   ! ==========================================================================
   ! Energy evaluation initial configuration
   ! ==========================================================================

   call recalculate_total_energy

   initial_energy_scaled = total_coulomb_energy*energy_conversion_factor*1.0d-3

   write (unit_output, '(a)') 'Initial configuration'
   write (unit_output, '(a, e12.5)') 'Coulomb energy (kJ/mol)  =', initial_energy_scaled

   ! ==========================================================================
   ! Start of Monte Carlo loop
   !
   ! Loop over mc_steps_inner*mc_steps_middle*mc_steps_outer particle moves.
   ! current_particle   = current particle
   ! ei    = new coulomb energy of particle il
   ! esa   = all coulomb interaction energies of the particle
   ! total_coulomb_energy = total energy (also ulaa(1))
   ! displacement_max(1..10)  displacement parameter for every ion
   !
   ! ==========================================================================

   do current_macro_step = 1, mc_steps_outer
      energy_accumulator_current(1) = 0.0

      do step_middle = 1, mc_steps_middle
         do step_inner = 1, mc_steps_inner
            current_particle = 1 + mod(total_mc_steps, num_particles)
            current_species = particle_species(current_particle)
            total_mc_steps = total_mc_steps + 1
            moves_per_species_total(current_species) = moves_per_species_total(current_species) + 1

            call random_number(random_value)
            trial_x = particle_x(current_particle) + displacement_max(current_particle)*(random_value - 0.5)
            call random_number(random_value)
            trial_y = particle_y(current_particle) + displacement_max(current_particle)*(random_value - 0.5)
            call random_number(random_value)
            trial_z = particle_z(current_particle) + displacement_max(current_particle)*(random_value - 0.5)

            if (trial_x > box_half) trial_x = trial_x - box_size
            if (trial_x < -box_half) trial_x = trial_x + box_size
            if (trial_y > box_half) trial_y = trial_y - box_size
            if (trial_y < -box_half) trial_y = trial_y + box_size
            if (trial_z > box_half) trial_z = trial_z - box_size
            if (trial_z < -box_half) trial_z = trial_z + box_size

            overlap_status = 99
            call evaluate_trial_move

            if (overlap_status == 0) then
               ! No hard core overlap, evaluate energy
               energy_difference = 0.0
               do loop_j = 1, num_particles
                  energy_difference = energy_difference + trial_energy(loop_j) - energy_matrix(loop_j, current_particle)
               end do

               boltzmann_factor = beta_inverse_temp*(energy_difference*energy_conversion_factor)

               ! Metropolis acceptance criterion
               if (boltzmann_factor < -80.0 .or. boltzmann_factor > 0.0) then
                  ! Automatically accept (very favorable or downhill move)
                  moves_per_species_accepted(current_species) = moves_per_species_accepted(current_species) + 1

                  ! Update particle position
                  particle_x(current_particle) = trial_x
                  particle_y(current_particle) = trial_y
                  particle_z(current_particle) = trial_z

                  ! Update energy matrix
                  do loop_j = 1, num_particles
                     energy_matrix(current_particle, loop_j) = trial_energy(loop_j)
                     energy_matrix(loop_j, current_particle) = trial_energy(loop_j)
                  end do

                  total_coulomb_energy = total_coulomb_energy + energy_difference
               else
                  ! Probabilistic acceptance for uphill moves
                  call random_number(random_value)
                  if (exp(boltzmann_factor) >= random_value) then
                     ! Accept move
                     moves_per_species_accepted(current_species) = moves_per_species_accepted(current_species) + 1

                     ! Update particle position
                     particle_x(current_particle) = trial_x
                     particle_y(current_particle) = trial_y
                     particle_z(current_particle) = trial_z

                     ! Update energy matrix
                     do loop_j = 1, num_particles
                        energy_matrix(current_particle, loop_j) = trial_energy(loop_j)
                        energy_matrix(loop_j, current_particle) = trial_energy(loop_j)
                     end do

                     total_coulomb_energy = total_coulomb_energy + energy_difference
                  else
                     ! Reject move due to energy
                     moves_per_species_energy_rejected(current_species) = moves_per_species_energy_rejected(current_species) + 1
                  end if
               end if
            else
               ! Reject move due to hard core overlap
               moves_per_species_hardcore_rejected(current_species) = moves_per_species_hardcore_rejected(current_species) + 1
            end if

            energy_accumulator_current(1) = energy_accumulator_current(1) + total_coulomb_energy

            if (mod(step_inner, widom_interval) == 0) then
               call calculate_collision_pressure1
               call calculate_widom_insertion1
            end if

         end do

         ! Accumulate RDF data at the end of each middle loop
         call accumulate_rdf
      end do

      energy_consistency_check = total_coulomb_energy
      call recalculate_total_energy

    if (total_coulomb_energy /= 0) energy_consistency_check = (energy_consistency_check - total_coulomb_energy)/total_coulomb_energy

      write (unit_macro, *)
      write (unit_macro, '(/, "************************************************************", &
      &"************************", /)')
      write (unit_macro, '(/ " MACROSTEP  ", i10, /, /, "Checkparameters : ", e10.3)') &
         current_macro_step, energy_consistency_check

      if (abs(energy_consistency_check) > 0.001) stop

      call recalculate_total_energy

      if (widom_interval <= mc_steps_inner) then
         call calculate_collision_pressure2
         call calculate_widom_insertion2
      end if

      energy_accumulator_current(1) = energy_accumulator_current(1)/dble(mc_steps_inner*mc_steps_middle)
      energy_accumulator_total(1) = energy_accumulator_total(1) + energy_accumulator_current(1)
      energy_accumulator_current(1) = energy_accumulator_current(1)*energy_conversion_factor*1.0d-3
      energy_per_macrostep(current_macro_step) = energy_accumulator_current(1)

      ! Calculate total acceptance and rejection
      total_accepted_moves = 0
      total_hardcore_rejections = 0
      total_energy_rejections = 0

      do loop_i = 1, num_species
         total_accepted_moves = total_accepted_moves + moves_per_species_accepted(loop_i)
         total_hardcore_rejections = total_hardcore_rejections + moves_per_species_hardcore_rejected(loop_i)
         total_energy_rejections = total_energy_rejections + moves_per_species_energy_rejected(loop_i)
      end do

      write (unit_macro, '(/, " total number of configurations    =", i8, /, &
      &" accepted configurations           =", i8, /, &
      &" energy rejected configurations    =", i8, /, &
      &" hard core rejected configurations =", i8, /)') &
         total_mc_steps, total_accepted_moves, total_energy_rejections, total_hardcore_rejections
      write (unit_macro, '(a, e12.5)') 'Coulomb energy (kj/mol)  =', energy_accumulator_current(1)
   end do

   ! ==========================================================================
   ! End of Monte Carlo loop
   !
   ! Calculation of the average energies and the pressures
   !
   ! ==========================================================================

   macro_step_inverse = 1.0/float(mc_steps_outer)
   energy_accumulator_total(1) = macro_step_inverse*energy_accumulator_total(1)
   energy_accumulator_total(1) = energy_accumulator_total(1)*energy_conversion_factor*1.0d-3

   call calculate_statistics(energy_per_macrostep, energy_variance, energy_accumulator_total(1), mc_steps_outer)

   write (unit_output, '(/)')
   write (unit_output, '(/, "************************************************************", &
   &"************************", /)')
   write (unit_output, '(/, a, i3, a, /)') 'FINAL RESULTS AFTER ', mc_steps_outer, &
      ' MACROSTEPS'
   write (unit_output, '(/, " total number of configurations    =", i8, /, &
   &" accepted configurations           =", i8, /, &
   &" energy rejected configurations    =", i8, /, &
   &" hard core rejected configurations =", i8, /)') &
      total_mc_steps, total_accepted_moves, total_energy_rejections, total_hardcore_rejections
   write (unit_output, '(/, "Divided per species", /, "species      accepted   ", &
   &" energy rej.   hard core rej.", /, 10(i6, 3(6x, f8.4), /), /)') &
      (loop_l, dble(moves_per_species_accepted(loop_l))/moves_per_species_total(loop_l), &
      dble(moves_per_species_energy_rejected(loop_l))/moves_per_species_total(loop_l), &
      dble(moves_per_species_hardcore_rejected(loop_l))/moves_per_species_total(loop_l), loop_l=1, num_species)
   write (unit_output, '(a, 2e12.5)') 'Coulomb energy (kj/mol)  =', energy_accumulator_total(1), energy_variance
   write (unit_output, *)

   if (widom_interval <= mc_steps_inner) then
      call calculate_collision_pressure3
      call calculate_widom_insertion3
   end if

   rewind unit_conf
   write (unit_conf, '(5e16.8)') (particle_x(loop_l), particle_y(loop_l), particle_z(loop_l), loop_l=1, num_particles)

   close (unit_conf)

   ! ==========================================================================
   !
   ! Calculation of pressures from the densities
   !
   ! ==========================================================================

   pressure_ideal = 0.0
   do loop_k = 1, num_species
      pressure_ideal = pressure_ideal + species_concentration(loop_k)*1.0d27/avogadro_number
   end do

   pressure_excess_energy = energy_accumulator_total(1)*1.0d30/(3*gas_constant*temperature*box_size**3*avogadro_number)
   pressure_excess_variance = energy_variance*1.0d30/(3*gas_constant*temperature*box_size**3*avogadro_number)
   pressure_total = pressure_ideal + pressure_collision_average + pressure_excess_energy
   pressure_total_variance  = sqrt(pressure_collision_variance* pressure_collision_variance+ pressure_excess_variance * pressure_excess_variance)

   write (unit_output, '(/, a, /)') 'TOTAL BULK PRESSURE '
   write (unit_output, '(a, f12.4, f10.4)') 'Ideal pressure       ', pressure_ideal, 0.0
   write (unit_output, '(a, f12.4, f10.4)') 'Energy con. <E/3V>   ', pressure_excess_energy, pressure_excess_variance
   write (unit_output, '(a, f12.4, f10.4)') 'Collision press      ', pressure_collision_average, pressure_collision_variance
   write (unit_output, '(70("_"))')
   write (unit_output, '(a, f12.4, f10.4)') 'Bulk pressure        ', pressure_total, pressure_total_variance

   ! Finalize and write RDF data
   call finalize_and_write_rdf

   stop

end program bulk

! =============================================================================
! =============================================================================
!
! Subroutine: earth
!
! Calculates the average and deviation of nnn values stored in
! array a(1..25)
!
! =============================================================================
! =============================================================================

subroutine calculate_statistics(data_array, standard_deviation, mean_value, array_size)

   implicit none
   integer :: i, array_size
   double precision :: data_array(25), standard_deviation, mean_value
   double precision :: sum_squared_deviations, normalization_factor

   normalization_factor = 1.0/(array_size*(array_size - 1))
   sum_squared_deviations = 0.0
   mean_value = 0.0

   do i = 1, array_size
      mean_value = mean_value + data_array(i)
   end do

   mean_value = mean_value/array_size

   do i = 1, array_size
      sum_squared_deviations = sum_squared_deviations + (data_array(i) - mean_value)*(data_array(i) - mean_value)
   end do

   sum_squared_deviations = sqrt(normalization_factor*sum_squared_deviations)
   standard_deviation = sum_squared_deviations

   return

end subroutine calculate_statistics

! =============================================================================
! =============================================================================
!
! Subroutine: earth2
!
! Calculates the average and deviation of nnn values stored in
! matrix a(1..25,1..num) and stores them in arrays xb and bbbwww
!
! =============================================================================
! =============================================================================

subroutine calculate_statistics_per_species(data_matrix, standard_deviations, mean_values, num_samples, num_species_to_process)

   implicit none
   integer, parameter :: max_species = 10
   integer :: i, j, num_samples, num_species_to_process
   double precision :: data_matrix(25, max_species), mean_values(max_species), standard_deviations(max_species)
   double precision :: current_std_dev, normalization_factor

   normalization_factor = 1.0/(num_samples*(num_samples - 1))

   do i = 1, num_species_to_process
      current_std_dev = 0.0
      mean_values(i) = 0.0

      do j = 1, num_samples
         mean_values(i) = mean_values(i) + data_matrix(j, i)
      end do

      mean_values(i) = mean_values(i)/num_samples

      do j = 1, num_samples
         current_std_dev = current_std_dev + (data_matrix(j, i) - mean_values(i))* &
                           (data_matrix(j, i) - mean_values(i))
      end do

      current_std_dev = sqrt(normalization_factor*current_std_dev)
      standard_deviations(i) = current_std_dev
   end do

   return

end subroutine calculate_statistics_per_species

! =============================================================================
! =============================================================================
!
! Subroutine: slump
!
! Set up of the random initial configuration.
!
! =============================================================================
! =============================================================================

subroutine initialize_random_configuration

   implicit none
   double precision :: random_value
   include 'bulk_f90.inc'

   ! Local variables
   integer :: i, particles_placed, placement_attempts, max_placement_attempts
   double precision :: delta_x, delta_y, delta_z, dist_squared
   double precision :: trial_position_x, trial_position_y, trial_position_z

   ! Note: Variable names have been updated to be more descriptive
   ! See bulk_f90.inc for the complete variable declarations

   write (unit_output, *) 'Random configuration generated'

   placement_attempts = 0
   max_placement_attempts = 100000
   particles_placed = 0

   do while (particles_placed < num_particles)
      placement_attempts = placement_attempts + 1

      if (placement_attempts > max_placement_attempts) then
         write (unit_output, *) ' too dense system'
         stop
      end if

      call random_number(random_value)
      trial_position_x = (random_value - 0.5)*box_size
      call random_number(random_value)
      trial_position_y = (random_value - 0.5)*box_size
      call random_number(random_value)
      trial_position_z = (random_value - 0.5)*box_size
      current_species = particle_species(particles_placed + 1)

      ! Check for overlaps with already placed particles
      overlap_status = 0
      do i = 1, particles_placed
         delta_x = trial_position_x - particle_x(i)
         delta_y = trial_position_y - particle_y(i)
         delta_z = trial_position_z - particle_z(i)
         delta_x = delta_x - aint(delta_x*box_half_inverse)*box_size
         delta_y = delta_y - aint(delta_y*box_half_inverse)*box_size
         delta_z = delta_z - aint(delta_z*box_half_inverse)*box_size
         dist_squared = delta_x**2 + delta_y**2 + delta_z**2

         if (dist_squared < hard_core_distance_sq(i, current_species)) then
            overlap_status = 1
            exit  ! Exit overlap check loop early
         end if
      end do

      ! If no overlap found, accept this particle position
      if (overlap_status == 0) then
         particles_placed = particles_placed + 1
         particle_x(particles_placed) = trial_position_x
         particle_y(particles_placed) = trial_position_y
         particle_z(particles_placed) = trial_position_z
      end if
   end do

end subroutine initialize_random_configuration

! =============================================================================
! =============================================================================
!
! Subroutine: qin
!
! Essential part of the program. Checks for overlap and calculate
! the energy change after a MC step.
!
! =============================================================================
! =============================================================================

subroutine evaluate_trial_move

   implicit none
   include 'bulk_f90.inc'

   ! Local variables
   integer :: i
   double precision :: delta_x, delta_y, delta_z, current_particle_charge
   double precision :: distance_squared(max_ions)  ! Squared distances for overlap checking

   ! Note: Variable names have been updated to be more descriptive
   ! See bulk_f90.inc for the complete variable declarations

   ! Compute all distances (vectorizable with branchless periodic boundary conditions)
   do i = 1, num_particles
      delta_x = trial_x - particle_x(i)
      delta_y = trial_y - particle_y(i)
      delta_z = trial_z - particle_z(i)

      ! Branchless periodic boundary conditions using aint
      delta_x = delta_x - aint(delta_x*box_half_inverse)*box_size
      delta_y = delta_y - aint(delta_y*box_half_inverse)*box_size
      delta_z = delta_z - aint(delta_z*box_half_inverse)*box_size

      distance_squared(i) = delta_x*delta_x + delta_y*delta_y + delta_z*delta_z
   end do

   ! Cancel out the hard core overlap with itself
   distance_squared(current_particle) = 1000000

   ! Check for hard core overlaps (early exit if found)
   do i = 1, num_particles
      if (distance_squared(i) < hard_core_distance_sq(i, current_species)) return
   end do

   overlap_status = 0
   current_particle_charge = particle_charge(current_particle)

   do i = 1, num_particles
      trial_energy(i) = current_particle_charge*particle_charge(i)/sqrt(distance_squared(i))
   end do

   trial_energy(current_particle) = 0

   return

end subroutine evaluate_trial_move

! =============================================================================
! =============================================================================
!
! Subroutine: liv
!
! Gives the total coulombic energy and the force by recalculating
! every particleinteraction. Used to check QIN
!
! =============================================================================
! =============================================================================

subroutine recalculate_total_energy

   implicit none
   include 'bulk_f90.inc'

   ! Local variables
   integer :: i, j, next_particle_index
   double precision :: delta_x, delta_y, delta_z, pairwise_energy
   double precision :: distance_squared(max_ions)  ! Squared distances for energy calculation

   ! Note: Variable names have been updated to be more descriptive
   ! See bulk_f90.inc for the complete variable declarations

   total_coulomb_energy = 0.0

   do i = 1, num_particles - 1
      trial_x = particle_x(i)
      trial_y = particle_y(i)
      trial_z = particle_z(i)
      next_particle_index = i + 1

      do j = next_particle_index, num_particles
         delta_x = trial_x - particle_x(j)
         delta_y = trial_y - particle_y(j)
         delta_z = trial_z - particle_z(j)
         delta_x = delta_x - aint(delta_x*box_half_inverse)*box_size
         delta_y = delta_y - aint(delta_y*box_half_inverse)*box_size
         delta_z = delta_z - aint(delta_z*box_half_inverse)*box_size
         distance_squared(j) = delta_x*delta_x + delta_y*delta_y + delta_z*delta_z
      end do

      ! Separate computation from reduction to enable vectorization
      do j = next_particle_index, num_particles
         pairwise_energy = particle_charge(i)*particle_charge(j)/sqrt(distance_squared(j))
         energy_matrix(i, j) = pairwise_energy
      end do

      ! Accumulate energy (separate reduction)
      do j = next_particle_index, num_particles
         total_coulomb_energy = total_coulomb_energy + energy_matrix(i, j)
      end do
   end do

   do i = 1, num_particles - 1
      do j = i + 1, num_particles
         energy_matrix(j, i) = energy_matrix(i, j)
      end do
   end do

   do i = 1, num_particles
      energy_matrix(i, i) = 0.0
   end do

   return

end subroutine recalculate_total_energy

! =============================================================================
!
! Subroutine: collision
!
! This routine checkes if a particle of type i in one halfcel is
! able to collide with a patricle of type j in the other. The
! contribution to the pressure is calculation by intergration
! At a random position around the real particle of type i a ghost
! particle of type j is inserted and the potential with every
! other particle in the system is calculated. This gives the
! collision probability
!
! =============================================================================
! =============================================================================

subroutine calculate_collision_pressure

   implicit none
   double precision :: random_value
   include 'bulk_f90.inc'

   ! Local variables
   integer :: i, j, k
   integer :: species_i, species_k
   integer :: num_particles_in_species, random_particle_index
   integer :: particle_list_start, total_collision_tests
   double precision :: relative_error(max_species)
   double precision :: collision_samples(25, max_species, max_species)
   double precision :: collision_matrix(max_species, max_species)
   double precision :: delta_x, delta_y, delta_z
   double precision :: axial_displacement, contact_distance, angle, radial_distance
   double precision :: ghost_x, ghost_y, ghost_z
   double precision :: total_electrostatic_potential
   double precision :: collision_avg, collision_var
   double precision :: temp_array(25)  ! Temporary array for statistical calculations
   double precision :: distance_squared(max_ions)  ! Squared distances for collision calculation
   double precision :: distance_inverse(max_ions)  ! Inverse distances for electrostatics

   ! Note: Variable names have been updated to be more descriptive
   ! See bulk_f90.inc for the complete variable declarations

   do i = 1, max_species
      do k = 1, max_species
         collision_matrix(i, k) = 0.0

         do j = 1, mc_steps_outer
            collision_samples(j, i, k) = 0.0
         end do
      end do
   end do

   return

   ! -------------------------------------------------------------------------

   entry calculate_collision_pressure1

   do i = 1, num_widom_insertions
      particle_list_start = 1

      do species_i = 1, num_species
         num_particles_in_species = int(species_properties(species_i, 2))

         do species_k = 1, num_species
            call random_number(random_value)
            random_particle_index = int(random_value*num_particles_in_species) + particle_list_start
            contact_distance = species_properties(species_i, 3) + species_properties(species_k, 3)
            call random_number(random_value)
            axial_displacement = (2*random_value - 1)*contact_distance
            call random_number(random_value)
            angle = 2*pi*random_value
            ghost_z = particle_z(random_particle_index) + axial_displacement
            radial_distance = sqrt(contact_distance*contact_distance - axial_displacement*axial_displacement) + 0.00001
            ghost_x = particle_x(random_particle_index) + radial_distance*cos(angle)
            ghost_y = particle_y(random_particle_index) + radial_distance*sin(angle)
            ghost_x = ghost_x - aint(ghost_x*box_half_inverse)*box_size
            ghost_y = ghost_y - aint(ghost_y*box_half_inverse)*box_size
            ghost_z = ghost_z - aint(ghost_z*box_half_inverse)*box_size

            ! Check for hard core overlaps and compute distances
            overlap_status = 0
            do k = 1, num_particles
               delta_x = dabs(ghost_x - particle_x(k))
               delta_y = dabs(ghost_y - particle_y(k))
               delta_z = dabs(ghost_z - particle_z(k))

               if (delta_x > box_half) delta_x = delta_x - box_size
               if (delta_y > box_half) delta_y = delta_y - box_size
               if (delta_z > box_half) delta_z = delta_z - box_size

               distance_squared(k) = delta_x*delta_x + delta_y*delta_y + delta_z*delta_z

               if (distance_squared(k) < hard_core_distance_sq(k, species_k)) then
                  overlap_status = 1
                  exit  ! Exit early if overlap found
               end if
            end do

            ! Only compute electrostatic potential if no hard core overlap
            if (overlap_status == 0) then
               do k = 1, num_particles
                  distance_inverse(k) = 1.0/sqrt(distance_squared(k))
               end do

               total_electrostatic_potential = 0

               do k = 1, num_particles
                  total_electrostatic_potential = total_electrostatic_potential + distance_inverse(k)*particle_charge(k)
               end do

               collision_matrix(species_i, species_k) = exp(beta_inverse_temp*total_electrostatic_potential* &
                                 species_properties(species_k, 4)*energy_conversion_factor) + collision_matrix(species_i, species_k)
            end if
         end do

         particle_list_start = particle_list_start + num_particles_in_species
      end do
   end do

   return

   ! -------------------------------------------------------------------------

   entry calculate_collision_pressure2

   total_collision_tests = num_widom_insertions*int(mc_steps_inner/widom_interval)*mc_steps_middle

   write (unit_macro, '(/, /, a)') 'COLLISION PRESSURE MATRIX'
   write (unit_macro, '(/, a, i6)') 'Total collision trials per species ', total_collision_tests
   write (unit_macro, '("Species          ", 10(12x, i6), /)') (i, i=1, num_species)

   do k = 1, num_species
      do i = 1, num_species
         collision_samples(current_macro_step, i, k) = collision_matrix(i, k)/total_collision_tests
         collision_matrix(i, k) = 0
      end do

      write (unit_macro, '("        pressure    ", 10e18.6, /)') (collision_samples(current_macro_step, i, k), i=1, num_species)
   end do

   return

   ! -------------------------------------------------------------------------

   entry calculate_collision_pressure3

   write (unit_output, '(/, a, /)') 'COLLISION MATRIX  AND RELATIVE ERROR'
   write (unit_output, '("Species          ", 10(12x, i6), /)') (i, i=1, num_species)

   do i = 1, num_species
      do k = 1, num_species
         collision_avg = 0

         do j = 1, mc_steps_outer
            contact_correlation(j, i, k) = collision_samples(j, i, k)
            temp_array(j) = collision_samples(j, i, k)
         end do

         call calculate_statistics(temp_array, collision_var, collision_avg, mc_steps_outer)

         relative_error(k) = 0
         if (collision_avg /= 0) relative_error(k) = collision_var/collision_avg
         collision_matrix(i, k) = species_concentration(i)*collision_avg
      end do

      write (unit_output, '(i4, "        ", 20e12.5)') i, (collision_matrix(i, k), relative_error(k), k=1, num_species)
   end do

   write (unit_output, '(/)')

   return

end subroutine calculate_collision_pressure

! =============================================================================
! =============================================================================
!
! Subroutine: widom
!
! Puts a ghost particle anywhere in the cell and calculates the
! energy with the rest of the system. Averaging of the exponent
! results in the chemical potential. A rescaling of the ion-
! charges is neccessary to assure neutrality. The integration
! is performed in entry 2 with a simpson's rule.
! The total chemical potential is divided in to the ideal (chid),
! hard core (chhc) and the electric contribution (chel). For
! each species the chem.pot. is given cellaverage, at the wall,
! and in the midplan.
!
! =============================================================================
! =============================================================================

subroutine calculate_widom_insertion

   implicit none
   double precision :: random_value
   include 'bulk_f90.inc'

   ! Local variables
   ! max_widom_species accounts for species at multiple measurement locations (bulk, wall, midplane)
   integer, parameter :: max_widom_species = max_species*3  ! max_species * (max_measurement_locations + 1)
   integer :: i, j, k
   integer :: measurement_position, total_chemical_potentials, total_widom_tests
   integer :: num_particles_in_species, particle_list_start, species_measurement_index
   integer :: integration_point_index, integration_index, rejection_sum
   integer :: hardcore_rejection_count(max_widom_species), is_rejected(max_widom_species)
   integer :: widom_count_per_macrostep(25), hardcore_rejection_all(0:5)
   double precision :: chem_pot_electrostatic(25, max_widom_species)
   double precision :: chem_pot_hardcore(25, max_widom_species)
   double precision :: chem_pot_excess(25, max_widom_species)
   double precision :: chem_pot_total(25, max_widom_species)
   double precision :: chem_pot_excess_widom(25, max_widom_species)
   double precision :: chem_pot_diff_wall(25, max_widom_species)
   double precision :: chem_pot_diff_midplane(25, max_widom_species)
   double precision :: chem_pot_integration(max_widom_species, 11)
   double precision :: energy_weighted_numerator(max_widom_species, 11)
   double precision :: energy_weighted_denominator(max_widom_species, 11)
   double precision :: chem_pot_integration_avg(max_widom_species, 11)
   double precision :: exponential_widom_sum(max_widom_species)
   double precision :: chem_pot_ideal(max_widom_species)
   double precision :: chem_pot_elec_avg(max_widom_species), chem_pot_elec_var(max_widom_species)
   double precision :: chem_pot_hc_avg(max_widom_species), chem_pot_hc_var(max_widom_species)
   double precision :: chem_pot_ex_avg(max_widom_species), chem_pot_ex_var(max_widom_species)
   double precision :: chem_pot_tot_avg(max_widom_species), chem_pot_tot_var(max_widom_species)
   double precision :: chem_pot_widom_avg(max_widom_species), chem_pot_widom_var(max_widom_species)
   double precision :: delta_x, delta_y, delta_z
   double precision :: test_particle_x, test_particle_y, test_particle_z
   double precision :: energy_weighted, energy_weighted_exponential, energy_weighted_lambda
   double precision :: total_electrostatic_potential, total_inverse_distance
   double precision :: simpson_weight_1, simpson_weight_2, simpson_weight_4, distance_sum
   double precision :: pairwise_energy
   double precision :: correlation_avg, correlation_var
   double precision :: temp_array(25)  ! Temporary array for statistical calculations
   double precision :: distance_squared(max_ions)  ! Squared distances for Widom insertion
   double precision :: distance_inverse(max_ions)  ! Inverse distances for electrostatics
   character(len=80) :: location_description(max_species)

   ! Note: Variable names have been updated to be more descriptive
   ! See bulk_f90.inc for the complete variable declarations

   do measurement_position = 0, measurement_location
      do i = 1, num_species
         if (species_concentration(i) /= 0) then
            chem_pot_ideal(i + num_species*measurement_position) = &
               dlog(species_concentration(i)/avogadro_number*1.0d27)
         else
            chem_pot_ideal(i + num_species*measurement_position) = -77
         end if
      end do
   end do

   do i = 1, max_species
      do j = 1, 25
         chem_pot_electrostatic(j, i) = 0
         chem_pot_hardcore(j, i) = 0
         chem_pot_excess(j, i) = 0
         chem_pot_total(j, i) = 0
      end do

      chem_pot_elec_avg(i) = 0
      chem_pot_hc_avg(i) = 0
      chem_pot_ex_avg(i) = 0
      chem_pot_tot_avg(i) = 0
      exponential_widom_sum(i) = 0
      hardcore_rejection_count(i) = 0

      do j = 1, 11
         energy_weighted_denominator(i, j) = 0
         energy_weighted_numerator(i, j) = 0
         chem_pot_integration_avg(i, j) = 0
      end do
   end do

   do i = 1, 25
      widom_count_per_macrostep(i) = 0
   end do

   do i = 0, 5
      hardcore_rejection_all(i) = 0
   end do

   location_description(1) = 'Chemical potentials averaged over the whole cell'
   location_description(2) = 'Chemical potentials at the wall'
   location_description(3) = 'Chemical potentials at the midplane'

   return

   entry calculate_widom_insertion1

   widom_count_per_macrostep(current_macro_step) = widom_count_per_macrostep(current_macro_step) + 1

   do measurement_position = 0, measurement_location
      do k = 1, num_widom_insertions
         call random_number(random_value)
         test_particle_x = box_size*(random_value - 0.5)
         call random_number(random_value)
         test_particle_y = box_size*(random_value - 0.5)
         test_particle_z = 0

         if (measurement_position <= 1) then
            call random_number(random_value)
            test_particle_z = box_size*(random_value - 0.5)
            if (measurement_position == 1) test_particle_z = sign(box_half, test_particle_z)
         end if

         ! Compute distances with branchless periodic boundary conditions (vectorizable)
         do j = 1, num_particles
            delta_x = test_particle_x - particle_x(j)
            delta_y = test_particle_y - particle_y(j)
            delta_z = test_particle_z - particle_z(j)

            delta_x = delta_x - aint(delta_x*box_half_inverse)*box_size
            delta_y = delta_y - aint(delta_y*box_half_inverse)*box_size
            delta_z = delta_z - aint(delta_z*box_half_inverse)*box_size

            distance_squared(j) = delta_x*delta_x + delta_y*delta_y + delta_z*delta_z
         end do

         rejection_sum = 0

         do i = 1, num_species
            species_measurement_index = measurement_position*num_species + i
            is_rejected(species_measurement_index) = 0

            do j = 1, num_particles
               if (distance_squared(j) < hard_core_distance_sq(j, i)) &
                  is_rejected(species_measurement_index) = 1
            end do

            rejection_sum = rejection_sum + is_rejected(species_measurement_index)
         end do

         ! Skip this test particle if all species reject it
         if (rejection_sum == num_species) then
            hardcore_rejection_all(measurement_position) = hardcore_rejection_all(measurement_position) + 1
            cycle  ! Continue to next test particle insertion
         end if

         do j = 1, num_particles
            distance_inverse(j) = 1.0/sqrt(distance_squared(j))
         end do

         total_electrostatic_potential = 0
         total_inverse_distance = 0
         particle_list_start = 1

         do i = 1, num_species
            num_particles_in_species = int(species_properties(i, 2))
            distance_sum = 0

            do j = particle_list_start, particle_list_start + num_particles_in_species - 1
               distance_sum = distance_sum + distance_inverse(j)
            end do

            particle_list_start = particle_list_start + num_particles_in_species
            total_electrostatic_potential = total_electrostatic_potential + distance_sum*species_properties(i, 4)
            total_inverse_distance = total_inverse_distance + distance_sum
         end do

         total_electrostatic_potential = total_electrostatic_potential + pairwise_energy

         do i = 1, num_species
            species_measurement_index = measurement_position*num_species + i

            if (is_rejected(species_measurement_index) == 1) then
               ! Species rejected due to hard core overlap
               hardcore_rejection_count(species_measurement_index) = hardcore_rejection_count(species_measurement_index) + 1
            else
               ! Species accepted, compute Widom insertion chemical potential
               exponential_widom_sum(species_measurement_index) = exponential_widom_sum(species_measurement_index) + &
                              exp(beta_inverse_temp*total_electrostatic_potential*species_properties(i, 4)*energy_conversion_factor)

               do integration_point_index = 0, 10
                  integration_index = integration_point_index + 1
                  energy_weighted = species_properties(i, 4)* &
                                    (total_electrostatic_potential - integration_point_index*0.1*species_properties(i, 4)* &
                                     total_inverse_distance/num_particles)
                  energy_weighted_lambda = energy_weighted*integration_point_index*0.1
                  energy_weighted_exponential = exp(beta_inverse_temp*energy_conversion_factor*energy_weighted_lambda)
                  energy_weighted_denominator(species_measurement_index, integration_index) = &
                     energy_weighted_denominator(species_measurement_index, integration_index) + energy_weighted_exponential
                  energy_weighted_numerator(species_measurement_index, integration_index) = &
                     energy_weighted_numerator(species_measurement_index, integration_index) - &
                     energy_weighted*beta_inverse_temp*energy_conversion_factor*energy_weighted_exponential
               end do
            end if
         end do
      end do
   end do

   return

   entry calculate_widom_insertion2

   total_chemical_potentials = num_species*(measurement_location + 1)

   do i = 1, total_chemical_potentials
      do j = 1, 11
         if (energy_weighted_denominator(i, j) == 0) then
            write (unit_output, *) ' WIDOM DENOMINATOR EQUALS ZERO', i, j
         else
            chem_pot_integration(i, j) = &
               energy_weighted_numerator(i, j)/energy_weighted_denominator(i, j)
            energy_weighted_numerator(i, j) = 0
            energy_weighted_denominator(i, j) = 0
            chem_pot_integration_avg(i, j) = chem_pot_integration_avg(i, j) + &
                                             1.0/mc_steps_outer*chem_pot_integration(i, j)
         end if
      end do

      simpson_weight_4 = chem_pot_integration(i, 2) + chem_pot_integration(i, 4) + &
                         chem_pot_integration(i, 6) + chem_pot_integration(i, 8) + &
                         chem_pot_integration(i, 10)
      simpson_weight_2 = chem_pot_integration(i, 3) + chem_pot_integration(i, 5) + &
                         chem_pot_integration(i, 7) + chem_pot_integration(i, 9)
      simpson_weight_1 = chem_pot_integration(i, 1) + chem_pot_integration(i, 11)
      chem_pot_electrostatic(current_macro_step, i) = &
         1.0/30.0*(simpson_weight_1 + 2*simpson_weight_2 + 4*simpson_weight_4)
   end do

   total_widom_tests = widom_count_per_macrostep(current_macro_step)*num_widom_insertions
   total_widom_tests = num_widom_insertions*int(mc_steps_inner/widom_interval)*mc_steps_middle

   do i = 1, total_chemical_potentials
      hardcore_rejection_count(i) = hardcore_rejection_count(i) + &
                                    hardcore_rejection_all(int((i - 1)/num_species))
      chem_pot_hardcore(current_macro_step, i) = &
         -dlog(dble(total_widom_tests - hardcore_rejection_count(i))/total_widom_tests)
      chem_pot_excess_widom(current_macro_step, i) = &
         -dlog(exponential_widom_sum(i)/total_widom_tests)
      chem_pot_excess(current_macro_step, i) = &
         chem_pot_electrostatic(current_macro_step, i) + chem_pot_hardcore(current_macro_step, i)
      chem_pot_total(current_macro_step, i) = &
         chem_pot_excess(current_macro_step, i) + chem_pot_ideal(i)
      exponential_widom_sum(i) = 0
      hardcore_rejection_count(i) = 0
   end do

   do i = 1, num_species
      chem_pot_diff_wall(current_macro_step, i) = &
         chem_pot_total(current_macro_step, i) - chem_pot_excess(current_macro_step, num_species + i)
      chem_pot_diff_midplane(current_macro_step, i) = &
         chem_pot_total(current_macro_step, i) - chem_pot_excess(current_macro_step, 2*num_species + i)
   end do

   write (unit_macro, '(/, /, "CHEMICAL POTENTIALS IN UNITS OF KT, MEASURED", i8, " TIMES")') total_widom_tests

   do measurement_position = 0, measurement_location
      hardcore_rejection_all(measurement_position) = 0
      write (unit_macro, '(/, a)') location_description(measurement_position + 1)
      write (unit_macro, '(/, "Type  Ideal     Hard Core       El. Static       Excess", &
      &"         Total         Uncorrected Excess")')
      write (unit_macro, '(10(i3, f8.4, 5(f9.4, "       "), /))') &
         (i, chem_pot_ideal(i), chem_pot_hardcore(current_macro_step, i), &
          chem_pot_electrostatic(current_macro_step, i), chem_pot_excess(current_macro_step, i), &
          chem_pot_total(current_macro_step, i), chem_pot_excess_widom(current_macro_step, i), &
          i=measurement_position*num_species + 1, (measurement_position + 1)*num_species)
   end do

   return

   entry calculate_widom_insertion3

   total_widom_tests = 0

   do i = 1, mc_steps_outer
      total_widom_tests = total_widom_tests + widom_count_per_macrostep(i)*num_widom_insertions
   end do

   total_chemical_potentials = num_species*(measurement_location + 1)

   call calculate_statistics_per_species(chem_pot_electrostatic, chem_pot_elec_var, chem_pot_elec_avg, &
                                         mc_steps_outer, total_chemical_potentials)
   call calculate_statistics_per_species(chem_pot_excess_widom, chem_pot_widom_var, chem_pot_widom_avg, &
                                         mc_steps_outer, total_chemical_potentials)
   call calculate_statistics_per_species(chem_pot_hardcore, chem_pot_hc_var, chem_pot_hc_avg, &
                                         mc_steps_outer, total_chemical_potentials)
   call calculate_statistics_per_species(chem_pot_excess, chem_pot_ex_var, chem_pot_ex_avg, &
                                         mc_steps_outer, total_chemical_potentials)
   call calculate_statistics_per_species(chem_pot_total, chem_pot_tot_var, chem_pot_tot_avg, &
                                         mc_steps_outer, total_chemical_potentials)

   write (unit_output, '(a, /)') 'CONTACT CORRELATION g(r)'
   write (unit_output, '("Species      ", 10(i12, 12x))') (i, i=1, num_species)

   pressure_collision_average = 0.0
   pressure_collision_variance = 0.0

   do i = 1, num_species
      do k = 1, num_species
         do j = 1, mc_steps_outer
            temp_array(j) = contact_correlation(j, i, k)* &
                            exp(chem_pot_excess_widom(j, k))
         end do

         call calculate_statistics(temp_array, correlation_var, correlation_avg, mc_steps_outer)
         contact_correlation(11, i, k) = correlation_var
         contact_correlation(12, i, k) = correlation_avg
      end do

      write (unit_output, '(i4, "        ", 10(f12.5, f10.5))') i, &
         (contact_correlation(12, i, k), contact_correlation(11, i, k), &
          k=1, num_species)
   end do

   write (unit_output, '(/, a, /)') 'CONTACT PRESSURE MATRIX'
   write (unit_output, '("Species      ", 10(i12, 12x))') (i, i=1, num_species)

   do i = 1, num_species
      do k = 1, num_species
         do j = 1, mc_steps_outer
            temp_array(j) = 2.0/3.0*pi*species_concentration(i)* &
                            contact_correlation(j, i, k)* &
                            exp(chem_pot_excess_widom(j, k) + chem_pot_ideal(k))* &
                            (species_properties(i, 3) + species_properties(k, 3))**3
         end do

         call calculate_statistics(temp_array, correlation_var, correlation_avg, mc_steps_outer)
         contact_correlation(11, i, k) = correlation_var
         contact_correlation(12, i, k) = correlation_avg
         pressure_collision_average = pressure_collision_average + correlation_avg
         pressure_collision_variance = pressure_collision_variance + correlation_var*correlation_var
      end do

      write (unit_output, '(i4, "        ", 10(f12.5, f10.5))') i, &
         (contact_correlation(12, i, k), contact_correlation(11, i, k), &
          k=1, num_species)
   end do

   pressure_collision_variance = sqrt(pressure_collision_variance)

   write (unit_output, '(/, a, f12.5, f10.5, /)') 'Total Collision pressure   =        &
   &', pressure_collision_average, pressure_collision_variance

   write (unit_output, '(/, "************************************************************", &
   &"**********************************************")')
   write (unit_output, '(/, /, "INTEGRAND IN WIDOM CORRECTION FOR DIFFERENT CHARGE PARAMETERS")')
   write (unit_output, '("TYPE", 11("  ", f4.2, "  "), /)') ((i - 1.0)/10.0, i=1, 11)

   do i = 1, total_chemical_potentials
      write (unit_output, '(10(i3, 11f8.3), /)') i, &
         (chem_pot_integration_avg(i, j), j=1, 11)
   end do

   write (unit_output, '(/, /, "CHEMICAL POTENTIALS IN UNITS OF KT, MEASURED", i8, " TIMES")') total_widom_tests

   do measurement_position = 0, measurement_location
      write (unit_output, '(/, a)') location_description(measurement_position + 1)
      write (unit_output, '(/, "Type  Ideal     Hard Core       El. Static       Excess", &
      &"         Total         Uncorrected Excess")')
      write (unit_output, '(10(i3, f8.4, 5(f9.4, f7.4), /))') &
         (mod(i - 1, num_species) + 1, &
          chem_pot_ideal(i), &
          chem_pot_hc_avg(i), chem_pot_hc_var(i), &
          chem_pot_elec_avg(i), chem_pot_elec_var(i), &
          chem_pot_ex_avg(i), chem_pot_ex_var(i), &
          chem_pot_tot_avg(i), chem_pot_tot_var(i), &
          chem_pot_widom_avg(i), chem_pot_widom_var(i), &
          i=measurement_position*num_species + 1, num_species*(measurement_position + 1))
      write (unit_output, '("2:1 ", 10(f8.4), /)') (exp((2*chem_pot_ex_avg(i - 2) + &
                                                         chem_pot_ex_avg(i - 1))/3))
   end do

   return

end subroutine calculate_widom_insertion

! =============================================================================
! =============================================================================
!
! Subroutine: initialize_rdf
!
! Initializes the radial distribution function (RDF) arrays and parameters
!
! =============================================================================
! =============================================================================

subroutine initialize_rdf

   implicit none
   include 'bulk_f90.inc'

   ! Local variables
   integer :: i, j, k

   ! Set RDF parameters
   ! Maximum distance is half the box size (due to periodic boundaries)
   rdf_max_distance = box_half
   rdf_bin_width = rdf_max_distance / rdf_num_bins

   ! Initialize histogram to zero
   do i = 1, num_species
      do j = 1, num_species
         do k = 1, rdf_num_bins
            rdf_histogram(k, i, j) = 0.0
         end do
      end do
   end do

   ! Initialize sample counter
   rdf_samples = 0.0

   write (unit_output, '(a)') 'RDF analysis initialized'
   write (unit_output, '(a, i6)') '  Number of bins: ', rdf_num_bins
   write (unit_output, '(a, f10.4)') '  Bin width (Angstrom): ', rdf_bin_width
   write (unit_output, '(a, f10.4, /)') '  Maximum distance (Angstrom): ', rdf_max_distance

   return

end subroutine initialize_rdf

! =============================================================================
! =============================================================================
!
! Subroutine: accumulate_rdf
!
! Accumulates histogram data for the radial distribution function
! by calculating all pairwise distances in the current configuration
!
! =============================================================================
! =============================================================================

subroutine accumulate_rdf

   implicit none
   include 'bulk_f90.inc'

   ! Local variables
   integer :: i, j, bin_index, species_i, species_j
   double precision :: delta_x, delta_y, delta_z, distance

   rdf_samples = rdf_samples + 1.0

   ! Loop over all unique particle pairs
   do i = 1, num_particles - 1
      species_i = particle_species(i)

      do j = i + 1, num_particles
         species_j = particle_species(j)

         ! Calculate distance with periodic boundary conditions
         delta_x = particle_x(i) - particle_x(j)
         delta_y = particle_y(i) - particle_y(j)
         delta_z = particle_z(i) - particle_z(j)

         ! Apply minimum image convention
         delta_x = delta_x - aint(delta_x*box_half_inverse)*box_size
         delta_y = delta_y - aint(delta_y*box_half_inverse)*box_size
         delta_z = delta_z - aint(delta_z*box_half_inverse)*box_size

         distance = sqrt(delta_x*delta_x + delta_y*delta_y + delta_z*delta_z)

         ! Add to histogram if within maximum distance
         if (distance < rdf_max_distance) then
            bin_index = int(distance/rdf_bin_width) + 1

            if (bin_index <= rdf_num_bins) then
               ! Add to histograms
               ! For same species: count both directions (i->j and j->i) by adding 2.0
               ! For different species: add to both histogram(i,j) and histogram(j,i)
               if (species_i == species_j) then
                  rdf_histogram(bin_index, species_i, species_j) = &
                     rdf_histogram(bin_index, species_i, species_j) + 2.0
               else
                  rdf_histogram(bin_index, species_i, species_j) = &
                     rdf_histogram(bin_index, species_i, species_j) + 1.0
                  rdf_histogram(bin_index, species_j, species_i) = &
                     rdf_histogram(bin_index, species_j, species_i) + 1.0
               end if
            end if
         end if
      end do
   end do

   return

end subroutine accumulate_rdf

! =============================================================================
! =============================================================================
!
! Subroutine: finalize_and_write_rdf
!
! Normalizes the RDF histograms and writes them to a CSV file
!
! =============================================================================
! =============================================================================

subroutine finalize_and_write_rdf

   implicit none
   include 'bulk_f90.inc'

   ! Local variables
   integer :: i, j, k, unit_rdf
   integer :: num_particles_i, num_particles_j
   double precision :: r_lower, r_upper, r_mid, shell_volume, number_density_j
   double precision :: ideal_count, normalization_factor, g_r

   unit_rdf = 20

   if (rdf_samples < 1.0) then
      write (unit_output, '(a)') 'Warning: No RDF samples collected'
      return
   end if

   ! Open CSV file for writing
   open (unit_rdf, file='rdf.csv', form='formatted', status='replace')

   ! Write header - simple approach
   write (unit_rdf, '(a)', advance='no') 'r'
   do i = 1, num_species
      do j = 1, num_species
         write (unit_rdf, '(a, i0, a, i0, a)', advance='no') ',g_', i, '_', j, '(r)'
      end do
   end do
   write (unit_rdf, *)  ! End header line

   ! Calculate and write normalized RDF for each bin
   do k = 1, rdf_num_bins
      r_lower = (k - 1)*rdf_bin_width
      r_upper = k*rdf_bin_width
      r_mid = (r_lower + r_upper)/2.0

      ! Write distance in first column
      write (unit_rdf, '(es15.7)', advance='no') r_mid

      ! Calculate and write g(r) for each species pair
      do i = 1, num_species
         num_particles_i = int(species_properties(i, 2))

         do j = 1, num_species
            num_particles_j = int(species_properties(j, 2))

            ! Calculate shell volume
            shell_volume = 4.0*pi/3.0*(r_upper**3 - r_lower**3)

            ! Number density of species j
            ! For same species, exclude self-particle from density
            if (i == j) then
               number_density_j = (num_particles_j - 1.0)/(box_size**3)
            else
               number_density_j = num_particles_j/(box_size**3)
            end if

            ! Expected number of j particles in shell around an i particle
            ! (for an ideal gas)
            ideal_count = shell_volume*number_density_j

            ! Normalization factor
            normalization_factor = rdf_samples*num_particles_i*ideal_count

            ! Calculate normalized g(r)
            if (normalization_factor > 0.0) then
               g_r = rdf_histogram(k, i, j)/normalization_factor
            else
               g_r = 0.0
            end if

            write (unit_rdf, '(a, es15.7)', advance='no') ',', g_r
         end do
      end do

      write (unit_rdf, *)  ! End of data line
   end do

   close (unit_rdf)

   write (unit_output, '(/, a)') 'RDF analysis completed'
   write (unit_output, '(a, f12.0)') '  Number of configurations sampled: ', rdf_samples
   write (unit_output, '(a, /)') '  RDF data written to rdf.csv'

   return

end subroutine finalize_and_write_rdf

! =============================================================================
! Helper subroutine to initialize the random number generator
! This is needed because the common block contains a variable named random_seed
! which conflicts with the intrinsic random_seed subroutine name
! =============================================================================
subroutine initialize_rng(seed_value)
   implicit none
   integer, intent(in) :: seed_value
   integer :: seed_size
   integer, allocatable :: seed_array(:)

   ! Get the size of the seed array
   call random_seed(size=seed_size)

   ! Allocate and initialize the seed array
   allocate (seed_array(seed_size))
   seed_array = seed_value

   ! Set the seed
   call random_seed(put=seed_array)

   deallocate (seed_array)
end subroutine initialize_rng
